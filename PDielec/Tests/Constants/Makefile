PYTEST := pytest

.PHONY: all test verbose fast clean

all: test

test:
	$(PYTEST) .

verbose:
	$(PYTEST) -v .

# All Constants tests are fast (no I/O, no heavy computation)
fast:
	$(PYTEST) -v .

clean:
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null; true
	find . -name "*.pyc" -delete 2>/dev/null; true
