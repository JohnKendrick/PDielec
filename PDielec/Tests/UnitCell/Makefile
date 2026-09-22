PYTEST := pytest

.PHONY: all test verbose fast clean

all: test

test:
	$(PYTEST) .

verbose:
	$(PYTEST) -v .

# Fast tests: skip symmetry (spglib) and CIF I/O
fast:
	$(PYTEST) test_A_length_units.py test_B_lattice.py test_C_coordinates.py test_D_volume_density.py test_G_coordinate_io.py -v

clean:
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null; true
	find . -name "*.pyc" -delete 2>/dev/null; true
