[![Build Status](https://travis-ci.org/GantZA/LatentOrderBookModel.jl.svg?branch=master)](https://travis-ci.org/GantZA/LatentOrderBookModel.jl)
[![Coverage Status](https://coveralls.io/repos/github/GantZA/LatentOrderBookModel.jl/badge.svg?branch=master)](https://coveralls.io/github/GantZA/LatentOrderBookModel.jl?branch=master)

# A Model of the Latent-Order Book Implemented in Julia

## Authors
* Michael Gant
* Tim Gebbie

## Acknowledgements

We would like to thank Byron Jacobs and Chris Angstman for their assistance and advice with the formulation of a numerical solution of the SPDE necessary for our implementation. We would also like to thank Donovan Platt for his research into calibrating the latent order book model and the code patterns for the numerical solution and calibration techniques. Michael Gant would like to thank Tim Gebbie for his supervision which entailed numerous discussions, valuable advice and resources.


## Installation

```
julia> ] add https://github.com/GantZA/LatentOrderBookModel.jl
```

## Arguments

The `ReactionDiffusionPricePaths` struct is created by instantiating the struct with the relevant parameters. Once instantiated, the object can be called as a function, with a seed as the only argument and a price path will be generated.

There are 2 command line interfaces which can be accesed. The commands and arguments are detailed below.
1) Running `$ julia src/main.jl arg1 arg2 ...`
2) After building the executable, `lobm_exec arg1 arg2 ...`

Each interface uses the same positional arguments which are:
* SEED :: Integer - The seed used for any random number generation. This ensures that paths are reproducible. A value of -1 will generate and use a random seed.
* T :: Integer - Number of time periods that are simulated
* τ :: Integer - Number of time periods that the DTWR solver uses.
* initial_mid_price :: Float - The initial price that the simulation begins at
* n_spatial_points :: Integer - The number of discretized price points used to solve the SPDE
* boltz_const :: Float - The Boltzmann constant used in calculating the Boltzmann Potential
* sample_std :: Float - The sample standard deviation of the price path. Used to find an upper and lower bound for the discretized price grid.
* D :: Float - The Diffusion coefficient in the SPDE
* nu :: Float - The latent order cancellation rate
* σ :: Float - The scaling value for the Stochastic Drift term
* λ :: Float - Source Term function parameter 1
* μ :: Float - Source Term function parameter 2
## Example Usage

### Julia Terminal

```
julia> using LatentOrderBookModel
julia> rdpp = ReactionDiffusionPricePaths(T=2300, τ=10, initial_mid_price=238.745,
  n_spatial_points=501, β=2.0, sample_std=7.415, D=5.0, nu=0.001, σ=1.0, λ=1.0, μ=0.5)
julia> rdpp(seed=57)

```

### Shell
```
$ julia src/main.jl 45 2300 10 238.745 501 2 7.415 0.001 5.0 0.001 1.0 1.0 0.5
```

### Compile using PackageCompiler.jl (requires master version)

```
julia> ] dev PackageCompiler
julia> using PackageCompiler
julia> PackageCompiler.compile_package("ArgParser" ,"LatentOrderBookModel")

```
After the Julia image has been compiled
```
julia -J ~/.julia/dev/PackageCompiler/sysimg/sys.so src/main.jl 45 2300 10 238.745 501 2 7.415 0.001 5.0 0.001 1.0 1.0 0.5

```

### Build Executable using PackageCompiler.jl (requires master version)

```
$ cd src
$ julia
julia> using PackageCompiler
julia> build_executable("lobm_exec.jl", "lobm_exec")

```
After the executable has been successfully built
```
$ ./builddir/lobm_exec 45 2300 10 238.745 501 2 7.415 0.001 5.0 0.001 1.0 1.0 0.5
```
