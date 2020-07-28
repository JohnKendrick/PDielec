p2cif-test:	
	@echo "Testing p2cif ...."
	@bash command.sh ${PYTHON_EXE} ${PREADER_OPTIONS} > all.cif
	@diff all.ref.cif all.cif

p2cif-regenerate:
	@echo "Regenerating p2cif reference files...."
	@bash command.sh ${PYTHON_EXE} ${PREADER_OPTIONS} > all.ref.cif
