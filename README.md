## Acknowledgements

We would like to thank Byron Jacobs and Chris Angstman for their work in deriving the SPDE,  its numerical solution and related code patterns. We would also like to thank Donovan Platt for his research into calibrating the latent order book model and the code patterns for the numerical solution and calibration techniques. Michael Gant would like to thank Tim Gebbie for his supervision which entailed numerous discussions, valuable advice and resources.


## Installation

```
julia> ] add https://github.com/GantZA/LatentOrderBookModel.jl
```

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
julia -J ~/.julia/dev/PackageCompiler/sysimg/sys.so src/main.jl 45 2300 10 .001 5.0 0.001 1.0 0.5

```
