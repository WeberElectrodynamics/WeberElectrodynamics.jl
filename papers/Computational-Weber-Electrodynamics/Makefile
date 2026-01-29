.PHONY: all paper figures setup clean help preview preview-open

# Default: full reproducibility pipeline
all: setup figures preview paper

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

# Preview figure standalone (fast iteration, output to build/)
preview:
	cd figures && pdflatex -output-directory=../build preview.tex
	cp build/preview.pdf figures/

# Preview and open
preview-open: preview
	open figures/preview.pdf

# Clean build artifacts
clean:
	rm -rf build/*
	rm -f Computational-Weber-Electrodynamics.pdf
	rm -f figures/preview.pdf

help:
	@echo "make all          - Full build (setup + figures + preview + paper)"
	@echo "make setup        - Install Julia dependencies"
	@echo "make figures      - Generate TikZ figures from Julia"
	@echo "make paper        - Compile LaTeX paper"
	@echo "make preview      - Compile figure standalone (fast)"
	@echo "make preview-open - Compile and open figure preview"
	@echo "make clean        - Remove build artifacts"
