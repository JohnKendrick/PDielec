test:	
	@echo "Testing preader for Abinit...."
	@bash command.sh > command.csv
	@echo "Comparing csv file with reference...."
	@diff command.ref.csv command.csv
	@echo "OK"

regenerate:
	@echo "Regenerating reference files from Abinit...."
	@bash command.sh > command.ref.csv
