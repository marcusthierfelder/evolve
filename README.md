evolve
======

generic time integrator library for PDEs
till now there is euler, icn, rk4, rk45 (time step adjustment)

all you need is to set the grid size (1D,2D,3D,...), add variables
as you want into a list of variables and provide a function which
calculates the righthandsides and pass it to the library. Then 
you can integrate. An simple output routine for vtk is included.

there is also a vtk and a netcdf output routine ... very rudimental

example
=======
advection_1d:
    advection equation which outputs ascii files 
    (can be ploted with ygraph and ncview)

lorenz
    simple lorenz model (0d)