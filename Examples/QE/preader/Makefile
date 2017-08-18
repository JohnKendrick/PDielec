test:	
	@echo "Testing preader for QE...."
	@bash command.sh > command.csv
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files from QE...."
	@bash command.sh > command.ref.csv
