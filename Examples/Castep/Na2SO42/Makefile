test:	
	@echo "Calculating absorption for Na2(SO4)2 from CASTEP...."
	@bash command.sh > command.out
	@echo "Comparing output files with reference"
	@diff command.ref.out command.out
	@echo "Comparing csv files with reference"
	@diff command.ref.csv command.csv
	@echo "OK"

regenerate:
	@echo "Regenerating reference files for Na2(SO4)2 from CASTEP...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv

