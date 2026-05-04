PYTEST := pytest

.PHONY: all test verbose fast clean

all: test

test:
	$(PYTEST) .

verbose:
	$(PYTEST) -v .

# Fast tests: skip slow effective-medium iterations
fast:
	$(PYTEST) test_A_tensor_utilities.py test_B_geometry.py test_C_depolarisation.py test_D_permittivity_conversions.py test_E_euler_angles.py -v

clean:
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null; true
	find . -name "*.pyc" -delete 2>/dev/null; true
