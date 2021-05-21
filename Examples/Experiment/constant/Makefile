pdgui:	
	@echo "PDGui Calculating constant permittivity from experimental...."
	@../../../pdgui -nosplash experiment constant.exp results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating constant permittivity from experimental...."
	@../../../pdgui -nosplash experiment constant.exp results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing constant permittivity from experimental...."
	@../../../pdgui experiment constant.exp results.xlsx -script script.py

