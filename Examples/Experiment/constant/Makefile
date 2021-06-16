pdgui:	
	@echo "PDGui Calculating Single Crystal for constant experiment file...."
	@../../../pdgui -nosplash experiment constant.exp results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating Single Crystal for constant experiment file...."
	@../../../pdgui -nosplash experiment constant.exp results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing Single Crystal for constant experiment file...."
	@../../../pdgui vasp OUTCAR -script script.py
	@../../../pdgui experiment constant.exp -script script.py

