pdgui:	
	@echo "Calculating absorption for Na2(SO4)2 from PHONOPY...."
	@../../../pdgui -nosplash phonopy OUTCAR.born results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel phonopy OUTCAR.born
	@../../../checkexcel results.xlsx results.ref.xlsx
	@diff OUTCAR.nma OUTCAR.nma.ref

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from PHONOPY...."
	@../../../pdgui -nosplash phonopy OUTCAR.born results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel phonopy OUTCAR.born
	@mv results.xlsx results.ref.xlsx
	@mv OUTCAR.nma OUTCAR.nma.ref

pdgui-view:
	@echo "Regenerating reference files for Na2(SO4)2 from PHONOPY...."
	@../../../pdgui phonopy OUTCAR.born -script script.py

