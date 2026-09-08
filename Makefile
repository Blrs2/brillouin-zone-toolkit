.PHONY: install test reproduce lint clean

install:
	python -m pip install -e '.[plot]'

test:
	PYTHONPATH=src python -m unittest discover -s tests -v

reproduce:
	PYTHONPATH=src python scripts/reproduce.py

lint:
	ruff check src tests scripts

clean:
	find . -type d -name '__pycache__' -prune -exec rm -rf {} +
	find . -type d -name '*.egg-info' -prune -exec rm -rf {} +
