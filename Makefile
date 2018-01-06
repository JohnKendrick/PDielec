#
# Installation assumes there is a SCRIPTS variable
# Defined in the environment which contains the
# destination directory. 
# If necessary the following line can be uncommented and edited
# SCRIPTS = /home/software/Scripts
# By default if no SCRIPTS environment variable is present then the installation goes to ~/bin

SCRIPTS ?= ~/bin

default:	
		@echo "Too install pdielec and preader in the $(SCRIPTS) directory"
		@echo "Type 'make install'"
		@echo "Too perform all the tests in the Examples directory"
		@echo "Type 'make test'"
		@echo "For a subset of the tests"
		@echo "Type 'make test_pdielec'"
		@echo " or  'make test_preader'"

install:	
		cp pdielec $(SCRIPTS)
		cp preader $(SCRIPTS)
		mkdir -p $(SCRIPTS)/Python
		cp -r Python/*.py $(SCRIPTS)/Python

test:		test_preader test_pdielec

tests:		test_preader test_pdielec

test_preader:		
		@echo "Testing preader functionality....."
		@( cd Examples; make --no-print-directory test_preader )

test_pdielec:		
		@echo "Testing pdielec functionality (takes a while to run)."
		@( cd Examples; make --no-print-directory test_pdielec )

regenerate:		
		@echo "Regenerating all reference data"
		@( cd Examples; make --no-print-directory regenerate )

pylint:		
		@pylint pdielec Python/*.py

pylama:		
		@pylama -i E501,E221,C901 .

clean:		
		@( cd Examples; make --no-print-directory clean )

