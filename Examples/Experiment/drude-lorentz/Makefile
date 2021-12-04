pdgui:	
	@echo "PDGui Calculating Single Crystal for MgO experiment file...."
	@../../../pdgui -nosplash experiment mgo.exp results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

benchmark:	
	@echo "PDGui benchmark Single Crystal for MgO experiment file...."
	@../../../pdgui -nosplash experiment mgo.exp -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating Single Crystal for MgO experiment file...."
	@../../../pdgui -nosplash experiment mgo.exp results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing Single Crystal for MgO experiment file...."
	@../../../pdgui vasp OUTCAR -script script.py
	@../../../pdgui experiment mgo.exp -script script.py

