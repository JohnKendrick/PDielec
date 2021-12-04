#
# Installation assumes there is a SCRIPTS variable
# Defined in the environment which contains the
# destination directory. 
# If necessary the following line can be uncommented and edited
# SCRIPTS = /home/software/Scripts
# By default if no SCRIPTS environment variable is present then the installation goes to ~/bin

SCRIPTS ?= ~/bin

default:	
		@echo "To install pdgui and preader in the $(SCRIPTS) directory"
		@echo "Type 'make install'"
		@echo "To perform all the tests in the Examples directory"
		@echo "Type 'make test'"
		@echo "For a subset of the tests"
		@echo " or  'make test-preader'"
		@echo " or  'make test-pdgui'"
		@echo " or  'make test-p2cif'"
		@echo " or  'make test-vibanalysis'"
		@echo " or  'make test-cli'"
		@echo " or  'make regenerate'"
		@echo " or  'make regenerate-vibanalysis'"
		@echo " or  'make pypi'"
		@echo " or  'make pyinstaller'"
		@echo " or  'make clean'"

install:	
		cp -P preader    $(SCRIPTS)
		cp -P p1reader   $(SCRIPTS)
		cp -P pdgui      $(SCRIPTS)
		cp -P p2cif      $(SCRIPTS)
		cp -P pdcompare  $(SCRIPTS)
		cp -P graphdatagenerator  $(SCRIPTS)
		cp -P vibanalysis  $(SCRIPTS)
		mkdir -p $(SCRIPTS)/PDielec
		mkdir -p $(SCRIPTS)/PDielec/GUI
		cp -r PDielec/*.py $(SCRIPTS)/PDielec
		cp -r PDielec/GUI/*.py $(SCRIPTS)/PDielec/GUI/
		cp -r PDielec/GUI/*.png $(SCRIPTS)/PDielec/GUI/

test:		test-pdgui test-cli

tests:		test-pdgui test-cli

test-cli:	test-preader test-p2cif test-vibanalysis

tests-cli:	test-preader test-p2cif test-vibanalysis

regenerate:	regenerate-pdgui regenerate-vibanalysis

.PHONY:		test-pdgui
test-pdgui:		
		@echo "Testing pdgui functionality....."
		@( cd Examples; $(MAKE) --no-print-directory test-pdgui )

.PHONY:		benchmark
benchmark:		
		@echo "Benchmarking....."
		@( cd Examples; $(MAKE) --no-print-directory benchmark )

.PHONY:		test-p2cif
test-p2cif:		
		@echo "Testing preader functionality....."
		@( cd Examples; $(MAKE) --no-print-directory test-p2cif )

.PHONY:		test-preader
test-preader:		
		@echo "Testing preader functionality....."
		@( cd Examples; $(MAKE) --no-print-directory test-preader )

.PHONY:		test-vibanalysis
test-vibanalysis:		
		@echo "Testing vibanalysis functionality....."
		@( cd Examples; $(MAKE) --no-print-directory test-vibanalysis )

.PHONY:		regenerate-pdgui
regenerate-pdgui:		
		@echo "Regenerating all reference data for pdgui"
		@( cd Examples; $(MAKE) --no-print-directory pdgui-regenerate )

.PHONY:		regenerate-vibanalysis
regenerate-vibanalysis:		
		@echo "Regenerating all reference data for vibanalysis"
		@( cd Examples; $(MAKE) --no-print-directory vibanalysis-regenerate )

.PHONY:		clean
clean:		
		@echo "Cleaning distribution and tests"
		@( cd Examples; $(MAKE) --no-print-directory clean )

.PHONY:		pyinstaller
pyinstaller:
		pyinstaller pdgui.spec -y
		cp -r dist/pdgui/PyQt5/Qt/plugins/platforms dist/pdgui

.PHONY:		pypi
pypi:		
		@echo "Generating PyPi distribution"
		@( rm -rf build dist PDielec.egg-info; python setup.py sdist bdist_wheel )
