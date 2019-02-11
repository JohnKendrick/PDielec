test-preader:
	@$(MAKE) --no-print-directory -C Castep/preader preader-test
	@$(MAKE) --no-print-directory -C Vasp/preader preader-test
	@$(MAKE) --no-print-directory -C Phonopy/preader preader-test
	@$(MAKE) --no-print-directory -C Gulp/preader preader-test
	@$(MAKE) --no-print-directory -C Crystal/preader preader-test
	@$(MAKE) --no-print-directory -C AbInit/preader preader-test
	@$(MAKE) --no-print-directory -C QE/preader preader-test

test-pdgui:
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid pdgui
	@$(MAKE) --no-print-directory -C Castep/Isoleucine pdgui
	@$(MAKE) --no-print-directory -C Castep/MgO pdgui
	@$(MAKE) --no-print-directory -C Castep/Bubbles pdgui
	@$(MAKE) --no-print-directory -C Castep/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C Vasp/F-Apatite pdgui
	@$(MAKE) --no-print-directory -C Vasp/ZnO pdgui
	@$(MAKE) --no-print-directory -C Vasp/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C Gulp/calcite pdgui
	@$(MAKE) --no-print-directory -C Gulp/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C Crystal/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C Crystal/Leucine pdgui
	@$(MAKE) --no-print-directory -C Crystal/Quartz pdgui
	@$(MAKE) --no-print-directory -C Crystal/ZnO/Default pdgui
	@$(MAKE) --no-print-directory -C Crystal/ZnO/CPHF pdgui
	@$(MAKE) --no-print-directory -C Crystal/ZnO/NoEckart pdgui
	@$(MAKE) --no-print-directory -C AbInit/AlAs pdgui
	@$(MAKE) --no-print-directory -C AbInit/BaTiO3 pdgui
	@$(MAKE) --no-print-directory -C AbInit/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C QE/ZnO pdgui
	@$(MAKE) --no-print-directory -C QE/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C QE/Cocaine pdgui
	@$(MAKE) --no-print-directory -C Phonopy/ZnO pdgui
	@$(MAKE) --no-print-directory -C Phonopy/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C SizeEffects/BaTiO3 pdgui
	@$(MAKE) --no-print-directory -C SizeEffects/MgO pdgui
	@$(MAKE) --no-print-directory -C SizeEffects/ZnO pdgui
	@$(MAKE) --no-print-directory -C Mie/MgO pdgui
	@$(MAKE) --no-print-directory -C Mie/MgO_lognormal pdgui
	@$(MAKE) --no-print-directory -C ATR/AlAs pdgui
	@$(MAKE) --no-print-directory -C ATR/Na2SO42 pdgui

clean:
	@echo "Cleaning old results.xlsx files"
	@find . -name results.xlsx -exec rm -f {} \;

pdgui-regenerate:
	@$(MAKE) --no-print-directory -C Castep/Isoleucine pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/MgO pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/Bubbles pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C Vasp/F-Apatite pdgui-regenerate
	@$(MAKE) --no-print-directory -C Vasp/ZnO pdgui-regenerate
	@$(MAKE) --no-print-directory -C Vasp/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C Gulp/calcite pdgui-regenerate
	@$(MAKE) --no-print-directory -C Gulp/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C Crystal/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C Crystal/Leucine pdgui-regenerate
	@$(MAKE) --no-print-directory -C Crystal/Quartz pdgui-regenerate
	@$(MAKE) --no-print-directory -C Crystal/ZnO/Default pdgui-regenerate
	@$(MAKE) --no-print-directory -C Crystal/ZnO/CPHF pdgui-regenerate
	@$(MAKE) --no-print-directory -C Crystal/ZnO/NoEckart pdgui-regenerate
	@$(MAKE) --no-print-directory -C AbInit/AlAs pdgui-regenerate
	@$(MAKE) --no-print-directory -C AbInit/BaTiO3 pdgui-regenerate
	@$(MAKE) --no-print-directory -C AbInit/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C QE/ZnO pdgui-regenerate
	@$(MAKE) --no-print-directory -C QE/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C QE/Cocaine pdgui-regenerate
	@$(MAKE) --no-print-directory -C Phonopy/ZnO pdgui-regenerate
	@$(MAKE) --no-print-directory -C Phonopy/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C SizeEffects/BaTiO3 pdgui-regenerate
	@$(MAKE) --no-print-directory -C SizeEffects/MgO pdgui-regenerate
	@$(MAKE) --no-print-directory -C SizeEffects/ZnO pdgui-regenerate
	@$(MAKE) --no-print-directory -C Mie/MgO pdgui-regenerate
	@$(MAKE) --no-print-directory -C Mie/MgO_lognormal pdgui-regenerate
	@$(MAKE) --no-print-directory -C ATR/AlAs pdgui-regenerate
	@$(MAKE) --no-print-directory -C ATR/Na2SO42 pdgui-regenerate
