PYTHON ?= python3
RSCRIPT ?= Rscript

.PHONY: all r python test-r test-python compare test clean

all: r python

r:
	$(RSCRIPT) R/run.R

python:
	$(PYTHON) python/run.py

test-r:
	$(RSCRIPT) R/tests/run_tests.R

test-python:
	$(PYTHON) -m pytest -c python/pyproject.toml python/tests

compare:
	$(PYTHON) python/validate_cross_language.py

test: test-r test-python compare

clean:
	rm -rf results/R results/python
