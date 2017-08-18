test:	
	@echo "Calculating absorption for ZnO (default settings) from Crystal...."
	@bash command.sh > command.out
	@../../../../checkcsv command.ref.out command.out
	@../../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for ZnO (default settings) from Crystal ...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv
