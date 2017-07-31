test:	
	@echo "Calculating absorption for MgO from CASTEP...."
	@bash command.sh > command.out
	@echo "  Comparing output files with reference"
	@diff command.ref.out command.out
	@echo "  Comparing csv files with reference"
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for MgO from CASTEP...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv

