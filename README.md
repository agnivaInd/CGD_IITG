# CGD_IITG
Codes for Computational Gas Dynamics problems. These are the codes that I had developed as part of my assignments and course project in the Computational Gas Dynamics course (ME683) at IIT Guwahati, during the final semester of my BTech. Codes are written in a fairly modular manner. 

[Assignment 1] Folder 1D_Wave_Equation_Various_Schemes deals with the solution of 1D Wave Equation with periodic boundary conditions for different numerical schemes viz. FTFS, FTBS, Upwind, Lax-Friedrichs and Lax-Wendroff. Mesh consists of 1D equi-sized cells. Initial values are assigned to cells by averaging over the cell. FVM used.

[Assignment 2] Folder Multiple_Linear_Wave_Equation deals with the solution of multiple linear wave equations (m x m system of equations) and their solution using the Upwind scheme. The initial conditions are Riemann ICs. Cells are equi-sized. FVM used.

[Coding Quiz] Folder Exact_Godunov_Solver_Burgers_Equation deals with the solution of the Burgers equation, using the exact solution of the Riemann problem through the Godunov scheme. Cells are equi-sized. FVM used

[Assignment 3/Course Project] Folder EulerEq_RoeApprox_StegerWarming deals with solving the 1D Sod Shock tube problem, using three different schemes - Roe Approximate Without Entropy Fix, Roe Approximate With Entropy Fix and Steger Warming (Roe Schemes are based on linearization of the Jacobian Matrix, whereas, Steger Warming is a flux splitting scheme). Solves the Euler Equation. Cells are equi-spaced. Zero flux boundary condition implemented. FVM used.
