module mod_DDC_mpi
use mpi
use kdtree2_module
implicit none

! ==============================================================================
!                              GLOBAL VARIABLES
! ==============================================================================

integer,    parameter :: spkind = kind(1.0), dpkind = kind(1.0d0)
! choose single or double precision, below
integer,    parameter :: pDki = dpkind
real(pDki), parameter :: pi = 4.0_pDki * atan(1.0_pDki)
integer               :: dim = 1 ! NOTE: this is the hard-coded 1 spatial dimension

! ==============================================================================
!                                 MPI VARIABLES
! ==============================================================================

integer                             :: ierror, num_cores, my_rank, namelen
integer, parameter                  :: master = 0
integer, parameter                  :: tag = 8273
integer, dimension(mpi_status_size) :: mpi_status
character(mpi_max_processor_name)   :: procname

! these are custom MPI types
integer                             :: realType  ! this will correspond to pDki, above
integer                             :: mpi_pType ! this will correspond to the derived ParticleType, below

! ==============================================================================

! derived type for particles
! NOTE: this is hard-coded for 1D right now
type ParticleType
    real(pDki) :: loc ! real-valued spatial location
    real(pDki) :: mass ! real-valued mass
    logical    :: active ! whether a given particle is currently active
    logical    :: ghost(2) ! indicates whether a particle is a ghost particle to the left or right. i.e., (1, 0) => left, (0, 1) => right
    logical    :: jumped(2) ! indicates whether the particle jumped subdomains and which way it went. i.e., (1, 0) => left, (0, 1) => right
end type ParticleType

! sparse matrix type in coordinate format
type SparseMatType
    integer    :: row
    integer    :: col
    real(pDki) :: val
end type SparseMatType

! holds the results of the kD tree fixed radius search
type kDRS_ResultType
    integer                 :: num    ! number of nearest neighbors (size of idx and rad)
    integer,    allocatable :: idx(:) ! indices of nearest neighbors
    real(pDki), allocatable :: rad(:) ! distances to nearest neighbors
end type kDRS_ResultType

! all the subroutines are below
contains

! subroutine to initialize the random number generator seed from clock time
! source: https://gcc.gnu.org/onlinedocs/gcc-4.2.1/gfortran/RANDOM_005fSEED.html
subroutine init_random_seed()
    integer              :: i, n, clock
    integer, allocatable :: seed(:)

    call random_seed(size = n)
    allocate (seed(n))
    call system_clock(count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed

! function to generate a linearly-spaced array with n points ranging [low, high]
function linspace(low, high, n) result(vec)
    real(pDki), intent(in   ) :: low, high
    integer,    intent(in   ) :: n
    real(pDki)                :: vec(n)
    real(pDki)                :: stride
    integer                   :: i

    if (n < 2) then
        print *, 'ERROR: linspace requires n > 1'
        call exit
    else
        vec(1) = low
        vec(n) = high
        stride = (high - low) / real(n - 1, pDki)
        do i = 2, n - 1
            vec(i) = low + real(i - 1, pDki) * stride
        enddo
    endif
end function linspace

! create an MPI real type to correspond to the chosen real kind
subroutine create_MPIRealType(x, realType)
    real(pDki), intent(in   ) :: x
    integer,    intent(  out) :: realType

    ! ========================== LOCAL VARIABLES ===============================
    integer :: p, r ! precision and range
    integer :: ierror

    p = precision(x)
    r = range(x)

    CALL MPI_Type_create_f90_real(p, r, realType, ierror)
end subroutine create_MPIRealType

! create and commit the ParticleType as an MPI type
subroutine build_derived_pType(mesg_mpi_pType)
    integer, intent(  out) :: mesg_mpi_pType

    ! ========================== LOCAL VARIABLES ===============================
    integer, parameter        :: number = 5 ! number of fields in the ParticleType
    integer                   :: block_lengths(number)
    integer(mpi_address_kind) :: displacements(number)
    integer                   :: typelist(number)
    real(pDki)                :: r
    logical                   :: log
    integer                   :: errcode, ierr

    typelist = (/ realType, realType, mpi_logical, mpi_logical, mpi_logical /)
    block_lengths = (/ 1, 1, 1, 2, 2 /)
    displacements = (/ int(0, mpi_address_kind), int(sizeof(r), mpi_address_kind),&
                       int(2, mpi_address_kind) * sizeof(r),&
                       int(2, mpi_address_kind) * sizeof(r) +  sizeof(log),&
                       int(2, mpi_address_kind) * sizeof(r) + int(3, mpi_address_kind) * sizeof(log)/)

    call mpi_type_create_struct(number, block_lengths, displacements, &
                                typelist, mesg_mpi_pType, ierr)
    if (ierr /= 0 ) then
        print *, 'Error in type create: ', ierr
        call mpi_abort(mpi_comm_world, errcode, ierr)
    endif

    call mpi_type_commit(mesg_mpi_pType, ierr)
    if (ierr /= 0 ) then
        print *, 'Error in type commit: ', ierr
        call mpi_abort(mpi_comm_world, errcode, ierr)
    endif
end subroutine build_derived_pType

! moves particles via random walk diffusion
subroutine diffuse(np, D, dt, idxActive, p)
    integer,            intent(in   ) :: np ! number of particles
    real(pDki),         intent(in   ) :: D, dt ! diffusion coefficient and time step
    integer,            intent(in   ) :: idxActive(np) ! array containing the indices of the active particles
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    real(pDki)                        :: normvec(np) ! vector which will hold Normal(0, 1) values

    ! call N(0, 1) generator
    call box_mullerp(np, normvec)

    p(idxActive)%loc = p(idxActive)%loc + sqrt(2.0_pDki * D * dt) * normvec
end subroutine diffuse

! this subroutine uses the Box-Muller transform to generate N(0,1)
! random numbers from U(0,1) numbers
! https://goo.gl/DQgmMu
! Note: this polar formulation seems to be consistently ~20% faster than the
! version that uses trig functions
! source for polar version (and standard version):
! https://www.taygeta.com/random/gaussian.html
subroutine box_mullerp(n, z)
    integer,    intent(in   ) :: n ! size of random vector to be generated
    real(pDki), intent(  out) :: z(n)

    ! ========================== LOCAL VARIABLES ===============================
    integer                   :: j
    real(pDki)                :: w, x1, x2
    real(pDki)                :: rand(2)

    ! initialize the random seed, just in case
    call init_random_seed()

    do j = 1, n/2
        w = 1.0_pDki
        do while (w >= 1.0_pDki)
            call random_number(rand)
            x1 = 2.0_pDki * rand(1) - 1.0_pDki
            x2 = 2.0_pDki * rand(2) - 1.0_pDki
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0_pDki * log(w)) / w)
        z(2 * j - 1 : 2 * j) = (/x1 * w, x2 * w/)
    enddo

    if (mod(n, 2) /= 0) then
        w = 1.0_pDki
        do while (w >= 1.0_pDki)
            call random_number(rand)
            x1 = 2.0_pDki * rand(1) - 1.0_pDki
            x2 = 2.0_pDki * rand(2) - 1.0_pDki
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0_pDki * log(w)) / w)
        z(n) = x1 * w
    endif
end subroutine box_mullerp

! reflective lower boundary
subroutine reflectLow(low, idxActive, p)
    real(pDki),         intent(in   ) :: low ! lower spatial boundary
    integer,            intent(in   ) :: idxActive(:) ! indices of the active particles
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! if a particle has exited the boundary, bounce it back the distance
    ! by which it overstepped
    where (p(idxActive)%loc < low) p(idxActive)%loc  = 2.0_pDki * low - p(idxActive)%loc
end subroutine reflectLow

! reflective upper boundary
subroutine reflectHigh(high, idxActive, p)
    real(pDki),         intent(in   ) :: high ! upper spatial boundary
    integer,            intent(in   ) :: idxActive(:) ! indices of the active particles
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! if a particle has exited the boundary, bounce it back the distance
    ! by which it overstepped
    where (p(idxActive)%loc > high) p(idxActive)%loc = 2.0_pDki * high - p(idxActive)%loc
end subroutine reflectHigh

! this is the big communications step between cores/subdomains
! any particle that random-walked across subdomain boundaries is sent to the
! neighboring subdomain, and any particle that is within paddist of the boundary
! is also sent to the neighboring subdomain as a ghost particle
subroutine swapDomains(pVecSize, Nactive, idxActive, Dlims, paddist, p)
    integer,            intent(in   ) :: pVecSize ! total size of the local particle array
    integer,            intent(inout) :: Nactive ! number of active particles
    integer,            intent(inout) :: idxActive(:) ! indices of the active particles
    real(pDki),         intent(in   ) :: Dlims(2) ! subdomain boundary limits
    real(pDki),         intent(in   ) :: paddist ! ghost particle padding distance
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    integer, allocatable :: idx(:) ! these are all used for bookkeeping in sending and receiving (hopefully what they are should be apparent by usage)
    integer              :: active(Nactive)
    integer              :: idxBig(pVecSize)
    integer              :: numSendLeft, numSendRight, recvSizeRight, recvSizeLeft
    integer              :: nJumpedRight, nJumpedLeft, totRecv, nGhostRight, nGhostLeft
    integer              :: i
    type(ParticleType)   :: tmp_pLeft(Nactive), tmp_pRight(Nactive) ! particle arrays for sending
    type(ParticleType)   :: recv_pLeft(Nactive), recv_pRight(Nactive) ! particle arrays for receiving

    ! initialize all of the relevant working info
    p%jumped(1) = .false.
    p%jumped(2) = .false.
    p%ghost(1) = .false.
    p%ghost(2) = .false.
    numSendLeft = 0
    numSendRight = 0
    recvSizeRight = 0
    recvSizeLeft = 0
    nJumpedRight = 0
    nJumpedLeft = 0
    nGhostRight = 0
    nGhostLeft = 0
    idxBig = (/(i, i = 1, pVecSize)/)

    ! this array will be used as an argument in the particle array
    allocate(idx(Nactive))
    active = idxActive(1 : Nactive)

    ! NOTE: subdomain boundaries are [lower, upper)

    ! first send particles to the left
    if (my_rank > 0) then
        ! if it is outside of the lower boundary, tag it to be sent to the left
        where (p(active)%loc < Dlims(1))
            p(active)%jumped(1) = .true.
        endwhere
        nJumpedLeft = count(p(active)%jumped(1))
        idx(1 : nJumpedLeft) = pack(active, p(active)%jumped(1))
        ! put it in a temporary array for sending
        tmp_pLeft(1 : nJumpedLeft) = p(idx(1 : nJumpedLeft))

        ! tag the particles that will be sent as ghost particles
        where (p(active)%loc >= Dlims(1) .and. p(active)%loc < Dlims(1) + paddist)
            p(active)%ghost(1) = .true.
        endwhere
        nGhostLeft = count(p(active)%ghost(1))
        idx(1 : nGhostLeft) = pack(active, p(active)%ghost(1))
        numSendLeft = nJumpedLeft + nGhostLeft
        ! add the ghost particles to the temporary particles array
        tmp_pLeft(nJumpedLeft + 1 : numSendLeft) = p(idx(1 : nGhostLeft))

        ! turn off the ghost particle indicator for the particles that are staying
        p%ghost(1) = .false.
        ! deactivate the particles that left the domain
        where (p(active)%jumped(1)) p(active)%active = .false.

        ! send the number of particles to be sent to the left
        call mpi_send(numSendLeft, 1, mpi_integer, my_rank - 1, tag + my_rank - 1,&
                      mpi_comm_world, ierror)
    endif
    if (my_rank < num_cores - 1) then
        ! receive the number of particles to be received from the right
        call mpi_recv(recvSizeRight, 1, mpi_integer, my_rank + 1, tag + my_rank, &
                      mpi_comm_world, mpi_status, ierror)
    endif

    ! if there are > 0 particles to be sent, send them
    if (my_rank > 0 .and. numSendLeft > 0) then
        call mpi_send(tmp_pLeft(1 : numSendLeft), numSendLeft, mpi_pType, my_rank - 1,&
                      tag + my_rank - 1, mpi_comm_world, ierror)
    endif
    ! if there are > 0 particles to be received, receive them
    if (my_rank < num_cores - 1 .and. recvSizeRight > 0) then
        call mpi_recv(recv_pRight(1 : recvSizeRight), recvSizeRight, mpi_pType,&
                      my_rank + 1, tag + my_rank, mpi_comm_world, mpi_status, ierror)
    endif

    ! now send particles to the right (this works the same as above, and thus is
    ! not commented)
    if (my_rank < num_cores - 1) then
        where (p(active)%loc >= Dlims(2))
            p(active)%jumped(2) = .true.
        endwhere
        nJumpedRight = count(p(active)%jumped(2))
        idx(1 : nJumpedRight) = pack(active, p(active)%jumped(2))
        tmp_pRight(1 : nJumpedRight) = p(idx(1 : nJumpedRight))

        where (p(active)%loc < Dlims(2) .and. p(active)%loc >= Dlims(2) - paddist)
            p(active)%ghost(2) = .true.
        endwhere
        nGhostRight = count(p(active)%ghost(2))
        idx(1 : nGhostRight) = pack(active, p(active)%ghost(2))
        numSendRight = nJumpedRight + nGhostRight
        tmp_pRight(nJumpedRight + 1 : numSendRight) = p(idx(1 : nGhostRight))

        p%ghost(2) = .false.
        where (p(active)%jumped(2)) p(active)%active = .false.

        call mpi_send(numSendRight, 1, mpi_integer, my_rank + 1, tag + my_rank + 1,&
                      mpi_comm_world, ierror)
    endif
    if (my_rank > 0) then
        call mpi_recv(recvSizeLeft, 1, mpi_integer, my_rank - 1, tag + my_rank, &
                      mpi_comm_world, mpi_status, ierror)
    endif

    if (my_rank < num_cores - 1 .and. numSendRight > 0) then
        call mpi_send(tmp_pRight(1 : numSendRight), numSendRight, mpi_pType, my_rank + 1,&
                      tag + my_rank + 1, mpi_comm_world, ierror)
    endif
    if (my_rank > 0 .and. recvSizeLeft > 0) then
        call mpi_recv(recv_pLeft(1 : recvSizeLeft), recvSizeLeft, mpi_pType,&
                      my_rank - 1, tag + my_rank, mpi_comm_world, mpi_status, ierror)
    endif

    deallocate(idx)
    allocate(idx(pVecSize))

    ! idx will show where to fill the newly received particles into the non-active spaces in "p"
    idx = pack(idxBig, .not. p%active)

    totRecv = recvSizeLeft + recvSizeRight

    ! if > 0 particles were received, put them in the main particles array,
    ! first filling them into the empty spaces of the non-active particles
    if (totRecv > 0) then
        recv_pLeft(recvSizeLeft + 1 : totRecv) = recv_pRight(1 : recvSizeRight)
        do i = 1, totRecv
            p(idx(i)) = recv_pLeft(i)
            p(idx(i))%active = .true.
        enddo
    endif
    deallocate(idx)
end subroutine swapDomains

! write plotting information to file
subroutine write_plot(uname, filename, vec, initial, N, Nsteps)
    integer,           intent(in   ) :: uname ! unit number to write to
    character(*),      intent(in   ) :: filename ! filename to write to
    real(pDki),        intent(in   ) :: vec(:) ! the vector to be written (either mass or location)
    logical,           intent(in   ) :: initial ! true if writing the initial time step
    integer, optional, intent(in   ) :: N, Nsteps ! number of particles and number of time steps, for the header

    ! if it is the initial time step, write the header indicating the shape of the the array
    if (initial) then
        open (unit = uname, file = filename, action = 'write')
        write (uname, *) N, Nsteps
        close (unit = uname, status='keep')
    endif

    open (unit = uname, file = filename, action = 'write', status='old', access='append')
    write (uname, *) vec
    close (unit = uname, status='keep')
end subroutine write_plot

! write run time to file for the current loop realization
subroutine write_time(uname, filename, time, initial, Num_ens, Np_ens, Np_list)
    integer,           intent(in   ) :: uname ! unit number to write to
    character(*),      intent(in   ) :: filename ! filename to write to
    real(pDki),        intent(in   ) :: time ! run time for the current realization
    logical,           intent(in   ) :: initial ! true if writing the for the first time => write the header
    integer, optional, intent(in   ) :: Num_ens, Np_ens ! number of runs in the ensemble (for averaging) and number of members in the particle number ensemble
    integer, optional, intent(in   ) :: Np_list(:) ! the number of particles in each member of the particle number ensemble

    ! if it is the initial time step, write the header indicating the shape of the the array
    if (initial) then
        open (unit = uname, file = filename, action = 'write')
        write (uname, *) Num_ens, Np_ens
        write (uname, *) Np_list
        close (unit = uname, status='keep')
    endif

    open (unit = uname, file = filename, action = 'write', status='old', access='append')
    write (uname, *) time
    close (unit = uname, status='keep')
end subroutine write_time

! write error to file for the current loop realization
subroutine write_error(uname, filename, error, initial, Num_ens, Np_ens, Np_list)
    integer,           intent(in   ) :: uname ! unit number to write to
    character(*),      intent(in   ) :: filename ! filename to write to
    real(pDki),        intent(in   ) :: error ! error for the current realization
    logical,           intent(in   ) :: initial ! true if writing the for the first time => write the header
    integer, optional, intent(in   ) :: Num_ens, Np_ens ! number of runs in the ensemble (for averaging) and number of members in the particle number ensemble
    integer, optional, intent(in   ) :: Np_list(:) ! the number of particles in each member of the particle number ensemble

    ! if it is the initial time step, write the header indicating the shape of the the array
    if (initial) then
        open (unit = uname, file = filename, action = 'write')
        write (uname, *) Num_ens, Np_ens
        write (uname, *) Np_list
        close (unit = uname, status='keep')
    endif

    open (unit = uname, file = filename, action = 'write', status='old', access='append')
    write (uname, *) error
    close (unit = uname, status='keep')
end subroutine write_error

! solves y = A * x, when A is in sparse coordinate format
! source: https://www.it.uu.se/education/phd_studies/phd_courses/pasc/lecture-1
subroutine SP_matVecMult(A, x, n, y)
    type(SparseMatType), intent(in   ) :: A(:)
    real(pDki),          intent(in   ) :: x(:)
    integer,             intent(in   ) :: n ! number of entries in A
    real(pDki),          intent(  out) :: y(:)
    integer                            :: i

    y = 0.0_pDki

    do i = 1, n
        y(A(i)%row) = y(A(i)%row) + A(i)%val * x(A(i)%col)
    enddo
end subroutine SP_matVecMult

! this is the sparse matrix-based mass transfer algorithm (SPaM)
subroutine massTrans_kDMat(n, idxActive, cutdist, denom, p)
    integer,            intent(in   ) :: n ! number of particles
    integer,            intent(in   ) :: idxActive(:) ! indices of the active particles
    real(pDki),         intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    real(pDki),         intent(in   ) :: denom ! denominator of the exponential in the co-location probability density
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    type(SparseMatType), allocatable :: Emat(:) ! mass transfer matrix
    integer                          :: start(n), finish(n), Nclose ! used for building the mass transfer matrix
    real(pDki)                       :: tmpmass(n) ! temporary array for holding particle masses
    type(ParticleType)               :: tmp_p(n) ! temporary particle array for dealing with ghost particles
    integer                          :: idx(n) ! indexing array
    integer                          :: i
    integer                          :: nNotGhost, idxNotGhost(n), idxNotGhostTmp(n) ! indexing array for non-ghost particles
    logical                          :: logNotGhost(n) ! logical array for non-ghost particles

    tmp_p = p(idxActive(1 : n))

    ! build the pairwise distance matrix
    call build_DistmatSparse(n, cutdist, tmp_p, Emat, start, finish, Nclose)
    ! build the matrix of co-location probabilities
    call build_PmatSparse(n, denom, start, finish, Emat)
    ! build the mass transfer matrix
    call build_EmatSparse(n, Nclose, Emat)

    tmpmass = tmp_p%mass

    ! conduct the mass transfers with sparse matrix-vector multiplication
    call SP_matVecMult(Emat, tmpmass, Nclose, tmp_p%mass)

    ! only change the masses of non-ghost particles
    idx = (/(i, i = 1, n)/)
    logNotGhost = tmp_p%ghost(1) .or. tmp_p%ghost(2)
    logNotGhost = .not. logNotGhost
    nNotGhost = count(logNotGhost)
    idxNotGhostTmp(1 : nNotGhost) = pack(idx, logNotGhost)
    idxNotGhost(1 : nNotGhost) = pack(idxActive(1 : n), logNotGhost)

    p(idxNotGhost(1 : nNotGhost))%mass = tmp_p(idxNotGhostTmp(1 : nNotGhost))%mass

    deallocate(Emat)
end subroutine massTrans_kDMat

! this is the looping sparse mass transfer algorithm (PLSS)
! portions of this work the same as in massTrans_kDMat(), and thus are not commented
subroutine massTrans_kDLoop(n, idxActive, cutdist, factor, denom, p)
    integer,            intent(in   ) :: n ! number of particles
    integer,            intent(in   ) :: idxActive(:) ! indices of the active particles
    real(pDki),         intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    real(pDki),         intent(in   ) :: factor ! constant multiplying the co-location probability density
    real(pDki),         intent(in   ) :: denom ! denominator of the exponential in the co-location probability density
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    type(kDRS_ResultType), allocatable :: neighbors(:) ! holds the results of the kD tree fixed radius search
    integer                            :: i, j ! loop iterators
    integer                            :: length
    integer, allocatable               :: idx(:) ! temporary array for holding indices of nearby particles (according to kD search)
    real(pDki), allocatable            :: rad(:) ! temporary array for holding distances to nearby particles (according to kD search)
    real(pDki), allocatable            :: Ptot(:) ! used to normalize the probabilities to sum to 1
    real(pDki)                         :: rescale ! used to normalize the probabilities to sum to 1
    real(pDki)                         :: vs_ds ! normalized co-location probability for the i-j particle pair
    real(pDki)                         :: dm ! mass differential for the i-j particle pair
    integer                            :: numTrans ! these next few are used for indexing arguments
    integer                            :: jpart
    logical, allocatable               :: logTrans(:)
    integer, allocatable               :: idxBig(:), idxTrans(:) ! used as indexing arguments to the particle array
    type(ParticleType)                 :: tmp_p(n) ! see above for what everything below this point is
    integer                            :: idxtmp(n)
    integer                            :: nNotGhost, idxNotGhost(n), idxNotGhostTmp(n)
    logical                            :: logNotGhost(n)

    tmp_p = p(idxActive(1 : n))

    ! conduct the kD tree fixed-radius search
    call fixedRadSearch(n, cutdist, tmp_p, neighbors)

    ! the following vectors are made as large as they possibly need to be to avoid
    ! having to allocate/deallocate repeatedly
    length = maxval(neighbors%num)
    allocate(idx(length), rad(length), Ptot(length), logTrans(length),&
             idxBig(length), idxTrans(length))

    idxBig = (/(i, i = 1, length)/)
    idxTrans = 0

    ! nested loop for conducting the mass transfers
    do i = 1, n
        ! we are only doing mass transfers for particles with an index greater
        ! than i, in order to avoid doing them twice
        logTrans(1 : neighbors(i)%num) = neighbors(i)%idx > i
        numTrans = count(logTrans(1 : neighbors(i)%num))
        if (numTrans > 0) then
            idxTrans(1 : numTrans) = pack(idxBig(1 : neighbors(i)%num), logTrans(1 : neighbors(i)%num))

            ! get the relevant indices and radii
            idx(1 : numTrans) = neighbors(i)%idx(idxTrans(1 : numTrans))
            rad(1 : numTrans) = neighbors(i)%rad(idxTrans(1 : numTrans))

            ! normalize the probabilities to sum to 1
            Ptot(1 : neighbors(i)%num) = factor * exp(neighbors(i)%rad / denom)
            rescale = sum(Ptot(1 : neighbors(i)%num))
            Ptot(1 : numTrans) = Ptot(idxTrans(1 : numTrans))

            ! make the indexing array
            idxTrans(1 : numTrans) = idxBig(1 : numTrans)

            ! do the transfers in a randomized order to avoid numerical artifacts
            call shuffleInt(idxTrans(1 : numTrans))

            ! do the pairwise mass transfers
            do j = 1, numTrans
                jpart = idx(idxTrans(j))
                vs_ds = Ptot(idxTrans(j)) / rescale
                dm = 0.5_pDki * (tmp_p(i)%mass - tmp_p(jpart)%mass) * vs_ds
                ! don't allow negative masses
                tmp_p(i)%mass     = maxval((/     tmp_p(i)%mass - dm, 0.0_pDki /))
                tmp_p(jpart)%mass = maxval((/ tmp_p(jpart)%mass + dm, 0.0_pDki /))
            enddo
            idxTrans = 0
        endif
    enddo

    ! only change the masses of non-ghost particles
    idxtmp = (/(i, i = 1, n)/)
    logNotGhost = tmp_p%ghost(1) .or. tmp_p%ghost(2)
    logNotGhost = .not. logNotGhost
    nNotGhost = count(logNotGhost)
    idxNotGhostTmp(1 : nNotGhost) = pack(idxtmp(1 : n), logNotGhost)
    idxNotGhost(1 : nNotGhost) = pack(idxActive(1 : n), logNotGhost)

    p(idxNotGhost(1 : nNotGhost))%mass = tmp_p(idxNotGhostTmp(1 : nNotGhost))%mass

    deallocate(idx, rad, Ptot, logTrans, idxTrans)
end subroutine massTrans_kDLoop

! this subroutine builds the (squared) distance matrix
! NOTE: hard-coded for 1D
subroutine build_DistmatSparse(n, cutdist, p, Distmat, start, finish, Nclose)
    integer,                          intent(in   ) :: n ! number of particles
    real(pDki),                       intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    type(ParticleType),               intent(in   ) :: p(:) ! particle array
    type(SparseMatType), allocatable, intent(  out) :: Distmat(:) ! sparse distance matrix
    integer,                          intent(  out) :: start(n), finish(n) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat
    integer,                          intent(  out) :: Nclose ! total number of neighbor particles found by the kD search (i.e, the length of the vectors in Distmat)

    ! ========================== LOCAL VARIABLES ===============================
    type(kDRS_ResultType), allocatable :: neighbors(:) ! holds the results of the kD tree fixed radius search
    integer                            :: i ! loop iterator
    integer                            :: tmpstart ! temporary variable to handle the n+1^th calculation of start in the loop

    ! conduct the kD tree fixed-radius search
    ! NOTE: this returns squared distances between particles
    call fixedRadSearch(n, cutdist, p, neighbors)

    ! allocate Distmat to have length = total number of neighbors found by the kD search
    Nclose = sum(neighbors%num)
    allocate(Distmat(Nclose))

    ! fill in Distmat
    tmpstart = 1
    do i = 1, n
        start(i) = tmpstart
        finish(i) = start(i) - 1 + neighbors(i)%num
        Distmat(start(i) : finish(i))%col = i
        Distmat(start(i) : finish(i))%row = neighbors(i)%idx
        Distmat(start(i) : finish(i))%val = real(neighbors(i)%rad, pDki)
        tmpstart = finish(i) + 1
    enddo

    deallocate (neighbors)
end subroutine build_DistmatSparse

! perform the kD tree fixed-radius search
subroutine fixedRadSearch(n, cutdist, p, neighbors)
    integer,                            intent(in   ) :: n ! number of particles
    real(pDki),                         intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    type(ParticleType),                 intent(in   ) :: p(:) ! particle array
    type(kDRS_ResultType), allocatable, intent(  out) :: neighbors(:) ! results of the fixed-radius search

    ! ========================== LOCAL VARIABLES ===============================
    type(kdtree2), pointer :: tree ! this is the KD tree
    real(kdkind)           :: locs(n) ! locations array to be passed to kD tree
    real(kdkind)           :: r2 ! value of squared search radius for kD tree

    ! convert particle locations to kdkind
    locs = real(  p%loc, kdkind)
    ! squared search cutoff distance in kdkind
    r2   = real(cutdist, kdkind)**2

    ! build the KD tree and search it
    call maketree(tree, dim, n, locs)

    allocate (neighbors(n))

    ! this finds the closest mobile particles to each immobile particle
    ! NOTE: this search returns the SQUARED distance between two particles
    call FRsearch(n, tree, r2, neighbors)
    call kdtree2_destroy(tree)
end subroutine fixedRadSearch

! build the sparse matrix of co-location probabilities
subroutine build_PmatSparse(n, denom, start, finish, Pmat)
    integer,             intent(in   ) :: n ! number of entries in the matrix
    real(pDki),          intent(in   ) :: denom ! denominator of the exponential in the co-location probability density
    integer,             intent(in   ) :: start(n), finish(n) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat
    type(SparseMatType), intent(inout) :: Pmat(:) ! sparse matrix of co-location probabilities

    ! ========================== LOCAL VARIABLES ===============================
    integer :: i

    ! calculate the encounter probabilities
    ! NOTE: the leading constant term is neglected here, due to the normalization
    ! NOTE: Pmat%val is already squared separation distance because of the kD search
    Pmat%val = exp(Pmat%val / denom)

    ! normalize the columns of Pmat to sum to 1
    do i = 1, n
        Pmat(start(i) : finish(i))%val = Pmat(start(i) : finish(i))%val /&
                                         sum(Pmat(start(i) : finish(i))%val)
    enddo
end subroutine build_PmatSparse

! build the sparse mass-transfer matrix
! \vec Emat = \vec I - 0.5 * [diag(\vec Pmat \vec 1)  - \vec Pmat]
!           = \vec I - 0.5 * [diag(rowsum(\vec Pmat)) - \vec Pmat]
!           = \vec I - 0.5 * diag(rowsum(\vec Pmat)) + 0.5 * \vec Pmat
subroutine build_EmatSparse(n, Nclose, Emat)
    integer,             intent(in   ) :: n ! total number of particles
    integer,             intent(in   ) :: Nclose ! total number of neighbor particles found by the kD search (i.e, the length of the vectors in Distmat)
    type(SparseMatType), intent(inout) :: Emat(:) ! sparse mass-transfer matrix

    ! ========================== LOCAL VARIABLES ===============================
    integer    :: i
    integer    :: diag(n) ! linear indices of the diagonal elements of Pmat
    real(pDki) :: rowsum(n) ! array holding the row sums of Pmat

    ! compute the rowsums of Pmat
    rowsum = 0.0_pDki
    do i = 1, Nclose
        rowsum(Emat(i)%row) = rowsum(Emat(i)%row) + Emat(i)%val
        if (Emat(i)%row == Emat(i)%col) then
            diag(Emat(i)%row) = i
        endif
    enddo

    ! this is the I - 0.5 * diag(P * 1) step
    rowsum = 1.0_pDki - 0.5_pDki * rowsum

    ! this is the 0.5 * \vec Pmat step
    Emat%val = 0.5_pDki * Emat%val

    ! finally, add them together
    Emat(diag)%val =  Emat(diag)%val + rowsum
end subroutine build_EmatSparse

! randomly shuffles an array of integers
! source: https://www.rosettacode.org/wiki/Knuth_shuffle#Fortran
subroutine shuffleInt(a)
    integer, intent(inout) :: a(:)
    integer                :: i, randpos, temp
    real                   :: r

    do i = size(a), 2, -1
        call random_number(r)
        randpos = int(r * i) + 1
        temp = a(randpos)
        a(randpos) = a(i)
        a(i) = temp
    enddo
end subroutine shuffleInt

! this builds a KD tree
subroutine maketree(tree2, d, n, locs)
    type(kdtree2), pointer, intent(  out) :: tree2 ! this is the KD tree
    integer,                intent(in   ) :: d, n ! number of spatial dimensions, number of particles
    real(kdkind),           intent(in   ) :: locs(d, n) ! location array for particles, with dimension d x n (number of spatial dimensions x number of particles)

    ! build the tree
    tree2 => kdtree2_create(locs, dim = d, sort = .false., rearrange = .true.)
        ! currently don't see a need to sort, as false appears to be quicker
        ! rearrange = true also appears to be quicker
end subroutine maketree

! this searches an already built KD tree
subroutine FRsearch(n, tree, r2, neighbors)
    integer,                intent(in   ) :: n ! number of particles
    type(kdtree2), pointer, intent(in   ) :: tree ! the KD tree
    real(kdkind),           intent(in   ) :: r2 ! squared search radius
    type(kDRS_ResultType),  intent(  out) :: neighbors(:) ! holds the results of the kD tree fixed radius search

    ! ========================== LOCAL VARIABLES ===============================
    integer                           :: i
    type(kdtree2_result), allocatable :: results(:) ! results array from KD tree module
    integer                           :: num_alloc

    ! allocate results as big as it could possibly be
    ! there's probably a more memory-efficient way to do this
    num_alloc = n
    allocate (results(num_alloc))

    ! loop over all particles
    do i = 1, n
        ! the type of search used here finds all the particles within
        ! squared distance r2 from the i^th particle in the list
        ! the hard-coded 0 is the 'correlation time' of the search
        ! as far as i can tell, this means a particle, itself, is included the
        ! idx list, while 1 would leave the particle, itself, out of the list
        call kdtree2_r_nearest_around_point(tree, i, 0, r2, neighbors(i)%num, num_alloc, results)

        ! allocate these based on how many nearby particles were found
        allocate (neighbors(i)%idx(neighbors(i)%num), neighbors(i)%rad(neighbors(i)%num))

        ! put the results in the neighbors array
        neighbors(i)%idx = results(1 : neighbors(i)%num)%idx
        neighbors(i)%rad = results(1 : neighbors(i)%num)%dis
    enddo

    deallocate (results)
end subroutine FRsearch

end module mod_DDC_mpi
