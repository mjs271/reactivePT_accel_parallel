# Parallel (Fortran/MPI) Code

This directory contains the code necessary to run the test cases for the PLSS and SPaM algorithms.

A `makefile` and run script (`run.sh`) are included, but no guarantees are made that they will work properly for all system configurations. They should, however, inform the process of building and running the executable.

The current simulation parameters in `DDC_mpi.f90` are designed to run a simple simulation that should print information to the screen, and the run should complete in under one minute. This run will create the files `PDDC_times.txt` and `PDDC_Error.txt`.

The `plotting` directory contains two Matlab scripts which may be used for plotting error/run time information (`plot_time_error_MPI.m`, which reads the files `PDDC_times.txt` and `PDDC_Error.txt`) and for plotting an animation of particle positions and masses (`plot_DDC_MPI.m`). Plotting the mass vs. position animation will require uncommenting some blocks of code that are indicated in the `DDC_mpi.f90` file. If plotting is turned on in this way, the files `mass.txt` and `locs.txt` will be written to be read by the Matlab script.

The `kdtree2.f90` code is sourced from the [jmhodges/kdtree2](https://github.com/jmhodges/kdtree2.git) repository.
