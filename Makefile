#
# Installation assumes there is a SCRIPTS variable
# Defined in the environment which contains the
# destination directory. 
# If necessary the following line can be uncommented and edited
# SCRIPTS = /home/software/Scripts
# By default if no SCRIPTS environment variable is present then the installation goes to ~/bin

SCRIPTS ?= ~/bin

default:	
		@echo "Too install pdgui and preader in the $(SCRIPTS) directory"
		@echo "Type 'make install'"
		@echo "Too perform all the tests in the Examples directory"
		@echo "Type 'make test'"
		@echo "For a subset of the tests"
		@echo " or  'make test-preader'"
		@echo " or  'make test-pdgui'"
		@echo " or  'make pyinstaller'"

.PHONY:		install
install:	
		cp -P preader    $(SCRIPTS)
		cp -P pdgui      $(SCRIPTS)
		cp -P p2cif      $(SCRIPTS)
		cp -P pdcompare  $(SCRIPTS)
		cp -P graphdatagenerator  $(SCRIPTS)
		mkdir -p $(SCRIPTS)/PDielec
		mkdir -p $(SCRIPTS)/PDielec/GUI
		cp -r PDielec/*.py $(SCRIPTS)/PDielec
		cp -r PDielec/GUI/*.py $(SCRIPTS)/PDielec/GUI/
		cp -r PDielec/GUI/*.png $(SCRIPTS)/PDielec/GUI/

test:		test-pdgui test-cli

test-cli:	test-preader 

tests-cli:	test-preader 

.PHONY:		pdgui
test-pdgui:		
		@echo "Testing pdgui functionality....."
		@( cd Examples; make --no-print-directory test-pdgui )

.PHONY:		test-preader
test-preader:		
		@echo "Testing preader functionality....."
		@( cd Examples; make --no-print-directory test-preader )

.PHONY:		regenerate
regenerate:		
		@echo "Regenerating all reference data for pdgui"
		@( cd Examples; make --no-print-directory pdgui-regenerate )

.PHONY:		clean
clean:		
		@( cd Examples; make --no-print-directory clean )

.PHONY:		pyinstaller
pyinstaller:
		pyinstaller pdgui.spec -y
		cp -r dist/pdgui/PyQt5/Qt/plugins/platforms dist/pdgui
