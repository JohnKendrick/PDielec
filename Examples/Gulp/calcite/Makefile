test:	
	@echo "Calculating absorption for Calcite from GULP...."
	@bash command.sh > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Calcite from GULP...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv
