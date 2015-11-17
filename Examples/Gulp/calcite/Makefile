test:	
	@echo "Calculating absorption for Calcite from GULP...."
	@bash command.sh > command.out
	@echo "Comparing output csv file reference...."
	@diff command.ref.out command.out
	@echo "Comparing csv file with reference...."
	@diff command.ref.csv command.csv
	@echo "OK"

regenerate:
	@echo "Regenerating reference files for Calcite from GULP...."
	@bash command.sh > command.ref.out
	@mv command.csv command.ref.csv
