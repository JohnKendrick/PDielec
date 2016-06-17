test:
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid test
	@$(MAKE) --no-print-directory -C Castep/MgO test
	@$(MAKE) --no-print-directory -C Castep/Na2SO42 test
	@$(MAKE) --no-print-directory -C Vasp/F-Apatite test
	@$(MAKE) --no-print-directory -C Vasp/ZnO test
	@$(MAKE) --no-print-directory -C Vasp/Na2SO42 test
	@$(MAKE) --no-print-directory -C Gulp/calcite test
	@$(MAKE) --no-print-directory -C Crystal/Na2SO42 test
	@$(MAKE) --no-print-directory -C Crystal/Leucine test
	@$(MAKE) --no-print-directory -C Crystal/Quartz test
	@$(MAKE) --no-print-directory -C Crystal/ZnO/Default test
	@$(MAKE) --no-print-directory -C Crystal/ZnO/CPHF test
	@$(MAKE) --no-print-directory -C Crystal/ZnO/NoEckart test
	@$(MAKE) --no-print-directory -C AbInit/AlAs test
	@$(MAKE) --no-print-directory -C AbInit/BaTiO3 test
	@$(MAKE) --no-print-directory -C AbInit/Na2SO42 test
	@$(MAKE) --no-print-directory -C QE/ZnO test
	@$(MAKE) --no-print-directory -C QE/Na2SO42 test
	@$(MAKE) --no-print-directory -C QE/Cocaine test

clean:
	@echo "Cleaning old command.out files"
	@find . -name command.out -exec rm -f {} \;
	@echo "Cleaning old command.csv files"
	@find . -name command.csv -exec rm -f {} \;

regenerate:
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid regenerate
	@$(MAKE) --no-print-directory -C Castep/MgO regenerate
	@$(MAKE) --no-print-directory -C Castep/Na2SO42 regenerate
	@$(MAKE) --no-print-directory -C Vasp/F-Apatite regenerate
	@$(MAKE) --no-print-directory -C Vasp/ZnO regenerate
	@$(MAKE) --no-print-directory -C Vasp/Na2SO42 regenerate
	@$(MAKE) --no-print-directory -C Gulp/calcite regenerate
	@$(MAKE) --no-print-directory -C Crystal/Na2SO42 regenerate
	@$(MAKE) --no-print-directory -C Crystal/Leucine regenerate
	@$(MAKE) --no-print-directory -C Crystal/Quartz regenerate
	@$(MAKE) --no-print-directory -C Crystal/ZnO/Default regenerate
	@$(MAKE) --no-print-directory -C Crystal/ZnO/CPHF regenerate
	@$(MAKE) --no-print-directory -C Crystal/ZnO/NoEckart regenerate
	@$(MAKE) --no-print-directory -C AbInit/AlAs regenerate
	@$(MAKE) --no-print-directory -C AbInit/BaTiO3 regenerate
	@$(MAKE) --no-print-directory -C AbInit/Na2SO42 regenerate
	@$(MAKE) --no-print-directory -C QE/ZnO regenerate
	@$(MAKE) --no-print-directory -C QE/Na2SO42 regenerate
	@$(MAKE) --no-print-directory -C QE/Cocaine regenerate
