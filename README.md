# DAB_HE


A Julia implementation of the Double Absorbing Boundary (DAB) 
method applied in a high-order spectral element method for 
solving the Helmholtz equation in 2D. This implementation 
follows the methodology presented in [1]. The current version of 
the code relates to the setup presented in Section 6.1 of the 
paper. This code was developed during my time as a postdoctoral fellow in the Faculty of Aerospace Engineering at the Technion - Israel Institute of Technolgoy.

<br>

[1] S. Papadimitropoulos, D. Givoli, The Double Absorbing Boundary method for the Helmholtz equation, Applied Numerical Mathematics (2021).


## Overview


This project implements a numerical solver for the Helmholtz equation using a high-order spectral element method with the Double Absorbing Boundary method. The implementation is written in Julia and provides efficient numerical methods for solving wave propagation problems.


## Project Structure


The project consists of several Julia modules:

- `DAB_HE.jl`: Main solver implementation
- `Mesh2Drect.jl`: 2D rectangular mesh generation
- `basis.jl` and `basisx.jl`: Basis functions implementation
- `intNodes.jl` and `intWeights.jl`: Integration nodes and weights
- `femerror.jl`: Error computation
- Various integration modules (I1.jl through I6.jl, ID.jl)
- Supporting utilities for sorting, node checking, and function evaluation


## Dependencies


- Julia programming language
- SparseArrays
- LinearAlgebra


## Usage


1. Ensure you have Julia installed on your system
2. Install the required dependencies:
   ```julia
   using Pkg
   Pkg.add(["SparseArrays", "LinearAlgebra"])
   ```
3. Clone this repository
4. Run the main solver:
   ```julia
   include("DAB_HE.jl")
   ```
