test:	
	@echo "Calculating absorption for F-Apatite from VASP...."
	@bash command.sh > command.out
	@echo "Comparing output csv file reference...."
	@diff command.ref.out command.out
	@echo "Comparing csv file with reference...."
	@diff command.ref.csv command.csv
	@echo "OK"

regenerate:
	@echo "Regenerating reference files for F-Apataite from CASTEP...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv
