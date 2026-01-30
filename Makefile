# Top-level Makefile for WeberElectrodynamics

CWE_DIR = papers/Computational-Weber-Electrodynamics

.PHONY: cwe cwe-paper cwe-figures cwe-preview cwe-preview-open cwe-setup cwe-clean

# Full build
cwe:
	$(MAKE) -C $(CWE_DIR) all

cwe-paper:
	$(MAKE) -C $(CWE_DIR) paper

cwe-figures:
	$(MAKE) -C $(CWE_DIR) figures

cwe-preview:
	$(MAKE) -C $(CWE_DIR) preview

cwe-setup:
	$(MAKE) -C $(CWE_DIR) setup

cwe-clean:
	$(MAKE) -C $(CWE_DIR) clean
