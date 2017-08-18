test:	
	@echo "Testing preader for Phonopy...."
	@bash command.sh > command.csv
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files from Phonopy...."
	@bash command.sh > command.ref.csv
