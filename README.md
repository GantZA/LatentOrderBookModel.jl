[![Build Status](https://travis-ci.org/GantZA/LatentOrderBookModel.jl.svg?branch=master)](https://travis-ci.org/GantZA/LatentOrderBookModel.jl)
[![Coverage Status](https://coveralls.io/repos/github/GantZA/LatentOrderBookModel.jl/badge.svg?branch=master)](https://coveralls.io/github/GantZA/LatentOrderBookModel.jl?branch=master)

# A Model of the Latent-Order Book Implemented in Julia

## Authors
* Michael Gant
* Tim Gebbie

## Acknowledgements

We would like to thank Byron Jacobs and Chris Angstmann for their assistance and advice with the formulation of a numerical solution of the SPDE necessary for our implementation. We would also like to thank Donovan Platt for his research into calibrating the latent order book model and the code patterns for the numerical solution and calibration techniques. Michael Gant would like to thank Tim Gebbie for his supervision which entailed numerous discussions, valuable advice and resources.


## Installation

```
julia> ] add https://github.com/GantZA/LatentOrderBookModel.jl
```

## Arguments

The `ReactionDiffusionPricePaths` object is created by instantiating the struct with the relevant parameters. Once instantiated, the object can be called as a function, with a seed as the only argument and LOB densities, price paths, right-jump probabilities and left-jump probabilities will be generated.

There are 2 command line interfaces which can be accesed. The commands and arguments are detailed below.
1) Running `$ julia src/main.jl arg1 arg2 ...`
2) After building the executable, `lobm arg1 arg2 ...`

Each command line interface uses the same positional arguments which are:
* SEED :: Integer - The seed used for any random number generation. This ensures that paths are reproducible. A value of -1 will generate and use random seeds.
* num_paths :: Integer - The number of price paths to simulate
* T :: Integer - Number of time periods that are simulated
* p₀ :: Float - The initial mid-price that each simulation begins at
* M :: Integer - The number of discretized price points used to solve the PDE
* β :: Float - The Boltzmann constant used in calculating the Boltzmann Potential
* sample_std :: Float - The sample standard deviation of the price path. Used to find an upper and lower bound for the discretized price grid.
* D :: Float - The Diffusion coefficient in the PDE
* σ :: Float - The scaling value for the Stochastic Drift term
* nu :: Float - The latent order cancellation rate, should be set to 0.0 in the Simple LOB model 
* λ :: Float - Source Term function parameter 1
* μ :: Float - Source Term function parameter 2
## Example Usage

### Julia Terminal

```
julia> using LatentOrderBookModel
julia> rdpp = ReactionDiffusionPricePaths(num_paths=1 ,T=2300, p₀=238.745,
  M=500, β=2.0, sample_std=20.0, D=5.0, σ=4.0, nu=0.0, λ=1.0, μ=0.5)
julia> rdpp(45)

```

### Shell
```
$ julia src/main.jl 45 1 2300 238.745 500 2.0 20.0 5.0 4.0 0.0 1.0 0.5
```

### Compile using PackageCompiler.jl (requires master version)

```
julia> ] add PackageCompiler#master
$ make compile
```
After the Julia image has been compiled
```
julia -J ~/.julia/dev/PackageCompiler/sysimg/sys.so src/main.jl 45 1 2300 238.745 500 2.0 20.0 5.0 4.0 0.0 1.0 0.5

```

### Build Executable using PackageCompiler.jl (requires master version)

```
$ make build

```
After the executable has been successfully built
```
$ ./builddir/lobm 45 1 2300 238.745 500 2.0 20.0 5.0 4.0 0.0 1.0 0.5
```
