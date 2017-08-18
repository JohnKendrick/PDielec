test:	
	@echo "Calculating absorption for Na2(SO4)2 from Quantum Espresso...."
	@bash command.sh > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Na2(SO4)2 from Quantum Espresso...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv
