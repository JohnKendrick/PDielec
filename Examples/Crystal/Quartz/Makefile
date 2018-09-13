pdgui:	
	@echo "PDGui Calculating absorption for Quartz from Crystal...."
	@../../../pdgui -nosplash crystal qua_hf_2d_f_ir-int.out results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Quartz from Crystal...."
	@../../../pdgui -nosplash crystal qua_hf_2d_f_ir-int.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for Quartz from Crystal...."
	@../../../pdgui crystal qua_hf_2d_f_ir-int.out -script script.py

