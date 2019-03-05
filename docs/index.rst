.. 2D Navier Stokes Solver in Stream Function Vorticity Form documentation master file, created by
   sphinx-quickstart on Tue Mar  5 17:15:30 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to documentation for the 2D Navier Stokes Solver in Stream Function Vorticity Form
==========================================================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Features
--------

- Cartesian Vorticity form of the Equations
- Hemholtz Hodge Decompositin into a vorticity equation and a stream-function Poisson's equation.
- Approximate factorization with lagged variables for the vorticity equations
- Interior nodes discretized with central 2nd order accurate differences
- Tridiagonal vorticity updates in x and y for fast sweeps.
- Boundary terms handled by augmented (1 sided FD)
- Gauss Seidel wtih successive over-relaxation for the Poisson Update
  of the scalar stream function.
- Fortran 90 implementation


Build
------------

Build the executable:

- make currently requires gfortran
- assumes you have make tools.
- cd into the src directory for this project and run:

$    make



Run
------------

Run the included fifi.dat inputs example with

$  ./test fifi.dat test1.out

- ./test runs the executable
- fifi.dat selects the included input file_exists

  This file itself has some documentation.

- follow the prompts to select time step and SOR factor.


Results
------------
- plotting of primitive variables u, v, stored in u.dat and v.dat
- plotting of stream functions stored in strm-func.dat
- plotting of vorticity stored in vorticity.dat
- requires Numpy and Matplotlib
- Run the python file to plot outputs:

$ python 2DCavityPlotTest.py


Notes
------------

This version tested on:

- OSX
- Linux
