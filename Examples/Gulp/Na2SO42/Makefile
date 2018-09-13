pdgui:	
	@echo "PDGui Calculating absorption for Na2(SO4)2 from GULP...."
	@../../../pdgui -nosplash gulp na2so42.gout results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from GULP...."
	@../../../pdgui -nosplash gulp na2so42.gout results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for Na2(SO4)2 from GULP...."
	@../../../pdgui gulp na2so42.gout -script script.py

