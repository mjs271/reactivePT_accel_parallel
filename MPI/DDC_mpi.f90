program DDC_mpi
use mod_DDC_mpi
implicit none

! ==============================================================================
!                           SIMULATION PARAMETERS
! ==============================================================================

real(pDki), parameter :: xlims(2)  = (/0.0_pDki, 150.0_pDki/) ! global domain limits
real(pDki), parameter :: midpt     = (xlims(2) + xlims(1)) / 2.0_pDki ! midpoint of the domain (used for the initial condition)
real(pDki), parameter :: D0        = 1e0_pDki ! coefficient of total diffusion in the system
real(pDki), parameter :: kappa_RW  = 0.5_pDki ! amount of diffusion to be simulated by random walks
real(pDki), parameter :: DRW       = kappa_RW * D0 ! diffusion coefficient for the random walks
real(pDki), parameter :: DMT       = (1.0_pDki - kappa_RW) * D0 ! diffusion coefficient for the mass transfers
real(pDki), parameter :: dt        = 1e-1_pDki ! time step
integer,    parameter :: Nsteps    = 1e2 ! number of time steps to take
real(pDki), parameter :: maxtime   = real(Nsteps, pDki) * dt ! simulation end time (used for RMSE calculation)
real(pDki), parameter :: cut_const = 3.0_pDki ! multiplier for cutdist
real(pDki), parameter :: pad_const = 3.0_pDki ! multiplier for paddist
real(pDki), parameter :: cutdist   = cut_const * sqrt(4.0_pDki * D0 * dt) ! number of standard deviations for kD search
real(pDki), parameter :: paddist   = pad_const * sqrt(4.0_pDki * D0 * dt) ! number of standard deviations for subdomain pad
integer,    parameter :: Num_ens   = 4 ! number of realizations in the ensemble
integer,    parameter :: Num_DDCMT = 2 ! corresponds to the 2 mass transfer schemes (i.e., SPaM = 1, PLSS = 2)


! ==============================================================================
!                           PARTICLE LOOP ARRAYS
! ==============================================================================

integer, parameter :: Np_ens          = 3 ! size of the particle number list (below)
integer, parameter :: Np_list(Np_ens) = (/ 5e2, 1e3, 2e3 /) ! list of particle numbers to loop over
! integer, parameter :: Np_list(Np_ens) = (/ 25e1, 5e2, 1e3, 5e3, 1e4, 25e3, 5e4, 1e5 /)  ! list of particle numbers to loop over

! ==============================================================================
!                             GENERAL VARIABLES
! ==============================================================================

integer                         :: i, Np_idx, ens, tstep, DDC_MTflag ! loop iterators
real(pDki),         allocatable :: a(:) ! CPU warmup array

type(ParticleType), allocatable :: p(:) ! ParticleType array (these are the particles)

integer                         :: Np ! total number of particles in the simulation
integer                         :: NDom ! number of subdomains (set below to be the number of cores in the MPI pool)
integer                         :: rem ! modulo remainder (used for load balancing, below)
integer                         :: coreNp, pVecSize ! number of particles initially on a core and the (larger) size the particle vector will be allocated to
real(pDki)                      :: Dlims(2) ! subdomain boundary limits

real(pDki)                      :: tic, toc ! timer variables for tracking run time
real(pDki)                      :: error, totError ! RMSE at final time (core local, and total)

real(pDki)                      :: ds, factor, denom ! these are used for co-location probability calculations
integer,            allocatable :: idx(:), idxActive(:) ! these are used to get the indices of active particles and idxActive is used as an argument to the particle array
integer                         :: Nactive ! number of active particles

! ==============================================================================
!                             READ/WRITE VARIABLES
! ==============================================================================

! unit number and filename for writing the run time arrays
integer,      parameter :: uTimeDDC       = 11
character(*), parameter :: timeDDCName    = 'PDDC_times.txt'

! unit number and filename for writing the error arrays
integer,      parameter :: uErrDDC        = 21
character(*), parameter :: ErrDDCName     = 'PDDC_Error.txt'

! ==============================================================================
! Plotting Variables (uncomment this block if you want to plot positions/masses)
! ==============================================================================

! ! particle array to hold all particles on master processor
! type(ParticleType), allocatable :: masterP(:)
! integer                         :: num, begin ! used for indexing in the masterP array

! ! unit numbers and filenames for writing the mass/location arrays
! integer,      parameter :: uLoc     = 31
! integer,      parameter :: uMass    = 32
! character(*), parameter :: locName  = 'locs.txt'
! character(*), parameter :: massName = 'mass.txt'

! ==============================================================================

! initialize the random seed right off the bat
call init_random_seed()

! initialize the mpi session
call mpi_init(ierror)
call mpi_comm_size(mpi_comm_world, num_cores, ierror)
call mpi_comm_rank(mpi_comm_world, my_rank, ierror)
call mpi_get_processor_name(procname, namelen, ierror)

! create an mpi real type to correspond to pDki
call create_MPIRealType(D0, realType)
! create and commit a derived mpi type to correspond to the ParticleType
call build_derived_pType(mpi_pType)

! print summary information to screen
if (my_rank == master) then
    write (*, *)
    print *, '==================== Run summary information ===================='
    write (*, "((i2, 1x), 'particle count(s),')") Np_ens
    write (*, "((i2, 1x), 'realization(s)')") Num_ens
    write (*, "((i2, 1x), 'DDC MT mode(s)')") Num_DDCMT
    print *, '================================================================='
    write (*, *)
endif

! ============================== Wake up the CPU ===============================
if (my_rank == master) then
    print *, 'Wait a moment, this will make the CPU work out a bit first...'
    write (*, *)
endif
do i = 1, 1200
    allocate(a(i))
    a = 1.0_pDki
    a = sqrt(a * real(i, pDki))
    deallocate(a)
enddo

if (my_rank == master) then
    print *, '--- Starting runs now ---'
    write (*, *)
endif

do Np_idx = 1, Np_ens ! particle number loop

    ! try and load balance the initial particles per core
    Np = Np_list(Np_idx)
    rem = mod(Np, num_cores)
    coreNp = Np / num_cores
    if (my_rank < rem) then
        coreNp = coreNp + 1
    endif

    ! ==========================================================================
    !         ***Uncomment this if you want to plot positions/masses***
    ! ==========================================================================
    ! allocate(masterP(Np))
    ! ==========================================================================

    ! set the local domain boundaries
    Dlims(1) = real(my_rank, pDki) * (xlims(2) - xlims(1)) / real(num_cores, pDki)
    Dlims(2) = real(my_rank + 1, pDki) * (xlims(2) - xlims(1)) / real(num_cores, pDki)

    ! allocate the local p vectors some factor larger, to account for ghost
    ! particles and more particles in one domain than the other after random walks start

    ! NOTE: this is currently an ad hoc solution and makes them unnecessarily large
        ! there is likely a better way to do this. however, everything in the
        ! Engdahl, et al., "Accelerating and Parallelizing... "
        ! paper will definitely run using this
    pVecSize = coreNp + coreNp
    pVecSize = ceiling(coreNp + (2.0_pDki * paddist) / (Dlims(2) - Dlims(1)) * 1.2_pDki * coreNp)
    allocate(p(pVecSize), idx(pVecSize), idxActive(pVecSize))

    ! master array of indices, to be used for logical indexing
    idx = (/(i, i = 1, pVecSize)/)

    ! number of subdomains is equal to the number of cores
    NDom = num_cores

    ! these are for the co-location probability calculations
        ! while ds and factor are not strictly necessary, they are included to
        ! match the "Accelerating and Parallelizing... " manuscript
    ds     = (xlims(2) - xlims(1)) / Np
    factor = ds / sqrt(8.0_pDki * pi * DMT * dt)
    denom  = -8.0_pDki * DMT * dt

    do DDC_MTflag = 1, Num_DDCMT ! mass-transfer mode loop (SPaM = 1, PLSS = 2)
        do ens = 1, Num_ens ! ensemble loop

            ! randomly (uniformly) assign the initial particle positions between
            ! the subdomain limits
            call random_number(p%loc)
            p(1 : coreNp)%loc = Dlims(1) + (Dlims(2) - Dlims(1)) * p(1 : coreNp)%loc

            ! initialize the particle type variables
            p%active = .false.
            p(1 : coreNp)%active = .true.
            p%mass = 0.0_pDki
            p(1 : coreNp)%mass = 0.0_pDki

            ! Heaviside initial condition
            where (p(1 : coreNp)%loc >= midpt) p(1 : coreNp)%mass = 1.0_pDki

            ! ==================================================================
            !  ***Uncomment this block if you want to plot positions/masses***
            ! ==================================================================

            ! ! send all the particles to master for plotting write-out
            ! masterP(1 : coreNp) = p(1 : coreNp)
            ! begin = coreNp + 1
            ! do i = 1, num_cores - 1

            !     if (my_rank == i) then
            !         call mpi_send(coreNp, 1, mpi_integer, master, tag, mpi_comm_world, ierror)
            !     endif
            !     if (my_rank == master) then
            !         call mpi_recv(num, 1, mpi_integer, i, tag, mpi_comm_world, mpi_status, ierror)
            !     endif
            !     if (my_rank == i) then
            !         call mpi_send(p(1 : coreNp), coreNp, mpi_pType, master, tag, mpi_comm_world, ierror)
            !     endif
            !     if (my_rank == master) then
            !         call mpi_recv(masterP(begin : begin + num - 1), num, mpi_pType,&
            !                       i, tag, mpi_comm_world, mpi_status, ierror)
            !     endif
            !     begin = begin + num

            ! enddo

            ! ! write out the plot information with initial == true,
            ! ! to write the header with the shape of the arrays
            ! if (my_rank == master) then
            !     call write_plot( uLoc,  locName,  masterP%loc,  .true., Np, Nsteps + 1)
            !     call write_plot(uMass, massName,  masterP%mass, .true., Np, Nsteps + 1)
            ! endif

            ! ==================================================================

            ! start the clock for recording run time
            tic = mpi_wtime()

            do tstep = 1, Nsteps ! time stepping loop

                ! get indices of the currently active particles
                idxActive = 0
                Nactive = count(p%active)
                idxActive(1 : Nactive) = pack(idx, p%active)

                ! random walk diffusion
                call diffuse(Nactive, DRW, dt, idxActive(1 : Nactive), p)

                ! reflecting boundary conditions
                if (my_rank == master) then
                    call reflectLow(xlims(1), idxActive(1 : Nactive), p)
                endif
                if (my_rank == num_cores - 1) then
                    call reflectHigh(xlims(2), idxActive(1 : Nactive), p)
                endif

                ! swap any particles that have crossed subdomain boundaries or
                ! are ghost particles to the neighboring subdomain/processor
                call swapDomains(pVecSize, Nactive, idxActive, Dlims, paddist, p)

                ! now that particles have been swapped between subdomains,
                ! recalculate the indices of currently active particles
                idxActive = 0
                Nactive = count(p%active)
                idxActive(1 : Nactive) = pack(idx, p%active)

                ! do the mass transfers
                select case (DDC_MTflag)
                    case (1) ! SPaM
                        call massTrans_kDMat(Nactive, idxActive, cutdist, denom, p)
                    case (2) ! PLSS
                        call massTrans_kDLoop(Nactive, idxActive, cutdist, factor, denom, p)
                    case default
                        print *, '**** ERROR: Invalid transfer mode ****'
                end select

                ! deactivate the ghost particles, now that mass transfers have occurred
                where (p%ghost(1)) p%active = .false.
                where (p%ghost(2)) p%active = .false.

                ! once again, recalculate the indices of currently active particles
                ! (not strictly necessary for the simulation, but required for
                ! the error calculations)
                idxActive = 0
                Nactive = count(p%active)
                idxActive(1 : Nactive) = pack(idx, p%active)

                ! ==============================================================
                !***Uncomment this block if you want to plot positions/masses***
                ! ==============================================================

                ! ! send all the particles to master for plotting write-out
                ! masterP(1 : Nactive) = p(idxActive(1 : Nactive))
                ! begin = Nactive + 1
                ! do i = 1, num_cores - 1

                !     if (my_rank == i) then
                !         call mpi_send(Nactive, 1, mpi_integer, master, tag, mpi_comm_world, ierror)
                !     endif
                !     if (my_rank == master) then
                !         call mpi_recv(num, 1, mpi_integer, i, tag, mpi_comm_world, mpi_status, ierror)
                !     endif
                !     if (my_rank == i) then
                !         call mpi_send(p(idxActive(1 : Nactive)), Nactive, mpi_pType, master, tag, mpi_comm_world, ierror)
                !     endif
                !     if (my_rank == master) then
                !         call mpi_recv(masterP(begin : begin + num - 1), num, mpi_pType,&
                !                       i, tag, mpi_comm_world, mpi_status, ierror)
                !     endif
                !     begin = begin + num

                ! enddo

                ! if (my_rank == master) then
                !     call write_plot( uLoc,  locName,  masterP%loc, .false.)
                !     call write_plot(uMass, massName, masterP%mass, .false.)
                ! endif

                ! ==============================================================

            enddo ! time stepping loop

            ! stock the clock for recording run time
            toc = mpi_wtime()

            ! write to screen some summary info and the run time
            if (my_rank == master) then
                select case (DDC_MTflag)
                    case (1) ! SPaM
                        write (*, "('SPaM', (i2, 1x), '(MT mode ', (i1, 1x), '): N = ', (i7, 1x),&
                               &' Ensemble number = ', (i2, 1x))") NDom, DDC_MTflag, Np, ens
                        write (*, *) 'run time = ', toc - tic
                    case (2) ! PLSS
                        write (*, "('PLSS', (i2, 1x), '(MT mode ', (i1, 1x), '): N = ', (i7, 1x),&
                               &' Ensemble number = ', (i2, 1x))") NDom, DDC_MTflag, Np, ens
                        write (*, *) 'run time = ', toc - tic
                    case default
                        print *, '**** ERROR: Invalid transfer mode ****'
                end select
            endif

            ! write the run time information to file
            if (my_rank == master) then
                if (ens == 1 .and. Np_idx == 1 .and. DDC_MTflag == 1) then
                    call write_time(uTimeDDC, timeDDCName, toc - tic,&
                                    .true., Num_ens, Np_ens, Np_list)
                else
                    call write_time(uTimeDDC, timeDDCName, toc - tic, .false.)
                endif
            endif

            ! compute error from analytic solution

            ! first compute individual "MSEs" on each core
            ! NOTE: these are divided by the TOTAL particle number so are not
                ! actually local mean-squared-errors
            error = sum((0.5_pDki * erfc(-((p(idxActive(1 : Nactive))%loc - midpt)&
                                           / sqrt(4.0_pDki * D0 * maxtime))) - p(idxActive(1 : Nactive))%mass)**2) / Np

            ! sum all the individual core's "MSEs"
            call mpi_reduce(error, totError, 1, realType, mpi_sum, master, mpi_comm_world, ierror)

            if (my_rank == master) then
                ! take the sqrt of the global MSE to get RMSE
                totError = sqrt(totError)
                ! print the error to screen
                print *,'totError = ', totError
                ! write the error information to file
                if (ens == 1 .and. Np_idx == 1 .and. DDC_MTflag == 1) then
                    call write_error(uErrDDC, ErrDDCName, totError,&
                                    .true., Num_ens, Np_ens, Np_list)
                else
                    call write_error(uErrDDC, ErrDDCName, totError, .false.)
                endif
            endif
        enddo ! ensemble loop
    enddo ! mass-transfer mode loop
    deallocate(p, idx, idxActive)
enddo ! particle number loop

call mpi_finalize(ierror)

end program DDC_mpi
