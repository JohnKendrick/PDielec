test:	
	@echo "Calculating absorption for ZnO (NOECKART settings) from Crystal...."
	@bash command.sh > command.out
	@echo "  Comparing output csv file reference...."
	@diff command.ref.out command.out
	@echo "  Comparing csv file with reference...."
	@../../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for ZnO (NOECKART settings) from Crystal ...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv
