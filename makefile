.PHONY: test
test:
	pytest src

.PHONY: setup
setup:
	conda create -n foundation-ingest python=3.6
	conda env update -n foundation-ingest -f environment.yml
