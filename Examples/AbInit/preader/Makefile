test:	
	@echo "Testing preader for Abinit...."
	@bash command.sh > command.csv
	@echo "  Comparing csv file with reference...."
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files from Abinit...."
	@bash command.sh > command.ref.csv
