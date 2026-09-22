PYTEST := pytest

.PHONY: all test verbose fast clean

all: test

test:
	$(PYTEST) .

verbose:
	$(PYTEST) -v .

# Fast tests: pure math and layer construction only
fast:
	$(PYTEST) test_A_math_utilities.py test_B_layer.py test_C_smatrix.py test_D_system.py -v

clean:
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null; true
	find . -name "*.pyc" -delete 2>/dev/null; true
