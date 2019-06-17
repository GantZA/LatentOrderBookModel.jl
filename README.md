[![Build Status](https://travis-ci.org/GantZA/LatentOrderBookModel.jl.svg?branch=master)](https://travis-ci.org/GantZA/LatentOrderBookModel.jl)

[![Coverage Status](https://coveralls.io/repos/github/GantZA/LatentOrderBookModel.jl/badge.svg?branch=master)](https://coveralls.io/github/GantZA/LatentOrderBookModel.jl?branch=master)

## Acknowledgements

We would like to thank Byron Jacobs and Chris Angstman for their work in deriving the SPDE,  its numerical solution and related code patterns. We would also like to thank Donovan Platt for his research into calibrating the latent order book model and the code patterns for the numerical solution and calibration techniques. Michael Gant would like to thank Tim Gebbie for his supervision which entailed numerous discussions, valuable advice and resources.


## Installation

```
julia> ] add https://github.com/GantZA/LatentOrderBookModel.jl
```

## Arguments

The `reaction_diffusion_path()` function can be evaluated in 3 ways:
1) directly in Julia
2) Running `src/main.jl`
3) executing `lobm_exec`.

Each method uses the same positional arguments which are:
* SEED :: Integer - The seed used for any random number generation. This ensures that paths are reproducible.
* T :: Integer - Number of time periods that are simulated
* τ :: Integer - Number of time periods that the DTWR solver uses.
* initial_mid_price :: Float - The initial price that the simulation begins at
* n_spatial_points :: Integer - The number of discretized price points used to solve the SPDE
* boltz_const :: Float - The Boltzmann constant used in calculating the Boltzmann Potential
* sample_std :: Float - The sample standard deviation of the price path. Used to find an upper and lower bound for the discretized price grid.
* σ :: Float - The standard deviation of the Advection Normal(mean=0) Random Variable.   
* D :: Float - The Diffusion coefficient in the SPDE
* η :: Float - The latent order cancellation rate.
* λ :: Float - Source Term function parameter 1
* μ :: Float - Source Term function parameter 2
## Example Usage

### Julia Terminal

```
julia> using LatentOrderBookModel
julia> reaction_diffusion_path(45, 2300, 10, 238.745, 501, 2, 7.415, 0.001, 5.0, 0.001, 1.0, 0.5)

```

### Shell
```
$ julia src/main.jl 45 2300 10 238.745 501 2 7.415 0.001 5.0 0.001 1.0 0.5
```

### Compile using PackageCompiler.jl (requires master version)

```
julia> ] dev PackageCompiler
julia> using PackageCompiler
julia> PackageCompiler.compile_package("ArgParser" ,"LatentOrderBookModel")

```
After the Julia image has been compiled
```
julia -J ~/.julia/dev/PackageCompiler/sysimg/sys.so src/main.jl 45 2300 10 238.745 501 2 7.415 0.001 5.0 0.001 1.0 0.5

```

### Build Executable using PackageCompiler.jl

```
$ cd src
$ julia
julia> using PackageCompiler
julia> build_executable("lobm_exec.jl", "lobm_exec")

```
After the executable has been successfully built
```
$ ./builddir/lobm_exec 45 2300 10 238.745 501 2 7.415 0.001 5.0 0.001 1.0 0.5
```
