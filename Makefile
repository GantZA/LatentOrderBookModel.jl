install:
	julia -e 'using Pkg; Pkg.add(PackageSpec(path="$(shell pwd)"))'
