test:
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid test
	@$(MAKE) --no-print-directory -C Castep/MgO test
	@$(MAKE) --no-print-directory -C Vasp/F-Apatite test
	@$(MAKE) --no-print-directory -C Vasp/ZnO test
	@$(MAKE) --no-print-directory -C Gulp/calcite test

clean:
	@echo "Cleaning old command.out files"
	@find . -name command.out -exec rm -f {} \;
	@echo "Cleaning old command.csv files"
	@find . -name command.csv -exec rm -f {} \;

regenerate:
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid regenerate
	@$(MAKE) --no-print-directory -C Castep/MgO regenerate
	@$(MAKE) --no-print-directory -C Vasp/F-Apatite regenerate
	@$(MAKE) --no-print-directory -C Vasp/ZnO regenerate
	@$(MAKE) --no-print-directory -C Gulp/calcite regenerate

