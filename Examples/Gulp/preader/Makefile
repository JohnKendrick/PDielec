test:	
	@echo "Testing preader for Gulp...."
	@bash command.sh > command.csv
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files from Gulp...."
	@bash command.sh > command.ref.csv
