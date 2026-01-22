.PHONY: all paper figures setup clean help

# Default: full reproducibility pipeline
all: setup figures paper

# Install Julia dependencies
setup:
	julia --project=code -e 'using Pkg; Pkg.instantiate()'

# Generate all figures
figures:
	julia --project=code code/example1_precession.jl

# Compile paper (output to build/, copy PDF to root)
paper:
	latexmk -pdf -output-directory=build Computational-Weber-Electrodynamics.tex
	cp build/Computational-Weber-Electrodynamics.pdf .

# Clean build artifacts
clean:
	rm -rf build/*
	rm -f Computational-Weber-Electrodynamics.pdf

help:
	@echo "make all     - Full build (setup + figures + paper)"
	@echo "make setup   - Install Julia dependencies"
	@echo "make figures - Generate TikZ figures from Julia"
	@echo "make paper   - Compile LaTeX paper"
	@echo "make clean   - Remove build artifacts"
