pdgui:	
	@echo "PDGui Calculating absorption for Calcite from GULP...."
	@../../../pdgui -nosplash gulp calcite.gout results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for Calcite from GULP...."
	@../../../pdgui -nosplash gulp calcite.gout -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Calcite from GULP...."
	@../../../pdgui -nosplash gulp calcite.gout results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for Calcite from GULP...."
	@../../../pdgui gulp calcite.gout -script script.py

