PYTEST := pytest

.PHONY: all test verbose fast clean

all: test

test:
	$(PYTEST) .

verbose:
	$(PYTEST) -v .

# Fast tests: skip database I/O (E4, E5 require spreadsheet)
fast:
	$(PYTEST) test_A_dielectric_constant.py test_B_dielectric_tabulated.py test_C_dielectric_lorentzian.py test_D_material_classes.py -v

clean:
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null; true
	find . -name "*.pyc" -delete 2>/dev/null; true
