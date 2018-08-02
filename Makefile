#
# Installation assumes there is a SCRIPTS variable
# Defined in the environment which contains the
# destination directory. 
# If necessary the following line can be uncommented and edited
# SCRIPTS = /home/software/Scripts
# By default if no SCRIPTS environment variable is present then the installation goes to ~/bin

SCRIPTS ?= ~/bin

default:	
		@echo "Too install pdielec, phonana and preader in the $(SCRIPTS) directory"
		@echo "Type 'make install'"
		@echo "Too perform all the tests in the Examples directory"
		@echo "Type 'make test'"
		@echo "For a subset of the tests"
		@echo "Type 'make test-pdielec'"
		@echo " or  'make test-preader'"
		@echo " or  'make test-phonana'"
		@echo " or  'make test-pdgui'"
		@echo " or  'make pyinstaller'"

.PHONY:		install
install:	
		cp pdielec.py $(SCRIPTS)
		cp preader.py $(SCRIPTS)
		cp phonana.py $(SCRIPTS)
		cp pdgui.py   $(SCRIPTS)
		cp p2cif.py   $(SCRIPTS)
		cp pdcompare.py  $(SCRIPTS)
		cp -P pdielec    $(SCRIPTS)
		cp -P preader    $(SCRIPTS)
		cp -P phonana    $(SCRIPTS)
		cp -P pdgui      $(SCRIPTS)
		cp -P p2cif      $(SCRIPTS)
		cp -P pdcompare  $(SCRIPTS)
		mkdir -p $(SCRIPTS)/Python
		mkdir -p $(SCRIPTS)/Python/PyMieScatt
		mkdir -p $(SCRIPTS)/Python/GUI
		cp -r Python/*.py $(SCRIPTS)/Python
		cp -r Python/PyMieScatt/*.py $(SCRIPTS)/Python/PyMieScatt/
		cp -r Python/GUI/*.py $(SCRIPTS)/Python/GUI/
		cp -r Python/GUI/*.png $(SCRIPTS)/Python/GUI/

test:		test-pdgui test-cli

test-cli:	test-preader test-pdielec test-phonana

tests-cli:	test-preader test-pdielec test-phonana

.PHONY:		pdgui
test-pdgui:		
		@echo "Testing pdgui functionality....."
		@( cd Examples; make --no-print-directory test-pdgui )

.PHONY:		test-phonana
test-phonana:		
		@echo "Testing phonana functionality....."
		@( cd Examples; make --no-print-directory test-phonana )

.PHONY:		test-preader
test-preader:		
		@echo "Testing preader functionality....."
		@( cd Examples; make --no-print-directory test-preader )

.PHONY:		test-pdielec
test-pdielec:		
		@echo "Testing pdielec functionality (takes a while to run)."
		@( cd Examples; make --no-print-directory test-pdielec )

.PHONY:		pdgui-regenerate
pdgui-regenerate:		
		@echo "Regenerating all reference data for pdgui"
		@( cd Examples; make --no-print-directory pdgui-regenerate )

.PHONY:		regenerate
regenerate:		
		@echo "Regenerating all reference data"
		@( cd Examples; make --no-print-directory regenerate )

.PHONY:		pylint
pylint:		
		@pylint pdielec Python/*.py

.PHONY:		pylama
pylama:		
		@pylama -i E501,E221,C901 .

.PHONY:		clean
clean:		
		@( cd Examples; make --no-print-directory clean )

.PHONY:		pyinstaller
pyinstaller:
		pyinstaller pdgui.spec -y
		cp -r dist/pdgui/PyQt5/Qt/plugins/platforms dist/pdgui
