test:	
	@echo "Calculating absorption for Na2(SO4)2 from VASP...."
	@bash command.sh > command.out
	@echo "  Comparing output csv file reference...."
	@diff command.ref.out command.out
	@echo "  Comparing csv file with reference...."
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Na2(SO4)2 from VASP...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv
