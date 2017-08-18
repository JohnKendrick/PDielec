test:	
	@echo "Testing preader for Vasp...."
	@bash command.sh > command.csv
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files from Vasp...."
	@bash command.sh > command.ref.csv
