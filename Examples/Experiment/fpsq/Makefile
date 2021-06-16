pdgui:	
	@echo "PDGui Calculating Single Crystal for fpsq experiment file...."
	@../../../pdgui -nosplash experiment fpsq.exp results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating Single Crystal for fpsq experiment file...."
	@../../../pdgui -nosplash experiment fpsq.exp results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing Single Crystal for fpsq experiment file...."
	@../../../pdgui vasp OUTCAR -script script.py
	@../../../pdgui experiment fpsq.exp -script script.py

