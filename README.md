# Differential-Equations
Some theory and numerical solvers for ODEs and PDEs
This repository contains several files regarding theoretical aspects and numerical solvers for differential equations. 
- FEM and Weak form of Elliptic Equations contains 2 folders. 1) Theory and Analysis which contains a PDF of examples proving the existence and uniqueness of elliptic PDEs and how weak formulations can be used for finite element method. 2) Finite Element Method contains a PDF with another example of deriving a linear system for FEM of a Sturm-Liouville differential equation. It also contains 2 folders: 1) with derivation and code for FEM solving the 1D heat equation, and 2) FEM method for 1D second order ode (with uniform and non-uniform mesh.
- Hyperbolic Problems contains 2 folders. 1) Contains code and derivation of the standard upwinding scheme for Burger's Equation. 2) Example of truncation error and stability analysis for hyperbolic PDEs. 
- Iterative Solvers contains code for Jacobi, Gauss-Siedel, and SOR (along with code for finding optimal parameter in SOR) solving a 2D boundary value problem, and a PDF with the results summarized.
- Nonlinear Boundary Value Problems contains 2 folders and a PDF. The PDF contains the derivation of Newton Iteration for two nonlinear boundary value problems and computational results. The two folders contain code for 1) nonlinear pendulum, and 2) a singular differential equation that exhibits boundary layer phenomena. 
- Asymptotics
- Finite Difference Solvers contains a folder with 2 PDFs and code. 1) PDF contains derivation of various finite differencing terms, and 2) PDF contains finite difference schemes for the heat equations -- both with centered differencing in space then 1) forward differencing in time, and 2) Crank-Nicholson in time.
- Additional Diff. Eqn. Theory contains a folder with PDF with examples of using the contraction mapping principle to establish the existence and uniqueness of solutions of differential equations and an example of the Implicit function theorem (for the same purpose). Green's function is under construction.

If I find the time, I would like to add some additional theory for linear PDEs such as transform theory, eigenfunction expansion, Sturm-Liouville Theory, and Green's Functions. 
