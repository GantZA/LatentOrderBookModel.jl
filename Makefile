install:
	julia -e 'using Pkg; Pkg.add(PackageSpec(path="$(shell pwd)"))'

compile:
	julia -e 'using PackageCompiler; PackageCompiler.compile_package("LatentOrderBookModel")'

build:
	julia -e 'using PackageCompiler; build_executable("lobm_exec.jl", "lobm")'
