ifndef VERBOSE
.SILENT:
endif

.PHONY: lint
lint:
	ruff check *.py tests/*.py
	mypy  *.py tests/*.py

.PHONY: format
format:
	ruff format *.py tests/*.py

.PHONY: test
test:
	pytest -v 