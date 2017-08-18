test:	
	@echo "Testing preader for Castep...."
	@bash command.sh > command.csv
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files from Castep...."
	@bash command.sh > command.ref.csv
