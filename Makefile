#
# Installation assumes there is a SCRIPTS variable
# Defined in the environment which contains the
# destination directory. 
# If necessary the following line can be uncommented and edited
# SCRIPTS = /home/software/Scripts


default:	
		@echo "Too install pdielec in the $(SCRIPTS) directory"
		@echo "Type 'make install'"
		@echo "Too perform all the tests in the Examples directory"
		@echo "Type 'make test'"

install:	
		cp pdielec $(SCRIPTS)
		mkdir -p $(SCRIPTS)/Python
		cp -r Python/*.py $(SCRIPTS)/Python

test:		
		@echo "The tests take a while to run."
		@echo "Verification is done by a simple diff with a reference file"
		@echo "Unfortunately just running on a different machine can cause small"
		@echo "and insignificant changes to the output files."
		@echo "If you see a change you must check if it is significant or not by hand."
		@( cd Examples; make --no-print-directory test )

clean:		
		@( cd Examples; make --no-print-directory clean )

