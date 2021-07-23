test-preader:
	@$(MAKE) --no-print-directory -C Castep/preader preader-test
	@$(MAKE) --no-print-directory -C Vasp/preader preader-test
	@$(MAKE) --no-print-directory -C Phonopy/preader preader-test
	@$(MAKE) --no-print-directory -C Gulp/preader preader-test
	@$(MAKE) --no-print-directory -C Crystal/preader preader-test
	@$(MAKE) --no-print-directory -C AbInit/preader preader-test
	@$(MAKE) --no-print-directory -C QE/preader preader-test

test-p2cif:
	@$(MAKE) --no-print-directory -C P2Cif/ p2cif-test

test-pdgui:
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid pdgui
	@$(MAKE) --no-print-directory -C Castep/Isoleucine pdgui
	@$(MAKE) --no-print-directory -C Castep/MgO pdgui
	@$(MAKE) --no-print-directory -C Castep/Bubbles pdgui
	@$(MAKE) --no-print-directory -C Castep/Na2SO42 pdgui
	@$(MAKE) --no-print-directory -C Castep/Castep17 pdgui
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
	@$(MAKE) --no-print-directory -C Experiment/constant pdgui
	@$(MAKE) --no-print-directory -C Experiment/fpsq pdgui
	@$(MAKE) --no-print-directory -C Experiment/drude-lorentz pdgui
	@$(MAKE) --no-print-directory -C Experiment/interpolation pdgui
	@$(MAKE) --no-print-directory -C Experiment/Mayerhofer pdgui
	@$(MAKE) --no-print-directory -C SingleCrystal/Bi2Se3 pdgui
	@$(MAKE) --no-print-directory -C SingleCrystal/Bi2Se3_film pdgui

test-vibanalysis:
	@$(MAKE) --no-print-directory -C VibAnalysis/AsparticAcid vibanalysis
	@$(MAKE) --no-print-directory -C VibAnalysis/BaTiO3 vibanalysis
	@$(MAKE) --no-print-directory -C VibAnalysis/Cocaine vibanalysis
	@$(MAKE) --no-print-directory -C VibAnalysis/F-Apatite vibanalysis
	@$(MAKE) --no-print-directory -C VibAnalysis/Isoleucine vibanalysis
	@$(MAKE) --no-print-directory -C VibAnalysis/Leucine vibanalysis

vibanalysis-regenerate:
	@$(MAKE) --no-print-directory -C VibAnalysis/AsparticAcid vibanalysis-regenerate
	@$(MAKE) --no-print-directory -C VibAnalysis/BaTiO3 vibanalysis-regenerate
	@$(MAKE) --no-print-directory -C VibAnalysis/Cocaine vibanalysis-regenerate
	@$(MAKE) --no-print-directory -C VibAnalysis/F-Apatite vibanalysis-regenerate
	@$(MAKE) --no-print-directory -C VibAnalysis/Isoleucine vibanalysis-regenerate
	@$(MAKE) --no-print-directory -C VibAnalysis/Leucine vibanalysis-regenerate

clean:
	@echo "Cleaning old results.xlsx files"
	@find . -name results.xlsx -exec rm -f {} \;
	@find . -name results.csv -exec rm -f {} \;
	@find . -name all.cif -exec rm -f {} \;
	@find . -name \*.nma -exec rm -f {} \;

pdgui-regenerate:
	@$(MAKE) --no-print-directory -C Castep/Isoleucine pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/AsparticAcid pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/MgO pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/Bubbles pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/Na2SO42 pdgui-regenerate
	@$(MAKE) --no-print-directory -C Castep/Castep17 pdgui-regenerate
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
	@$(MAKE) --no-print-directory -C Experiment/fpsq pdgui-regenerate
	@$(MAKE) --no-print-directory -C Experiment/constant pdgui-regenerate
	@$(MAKE) --no-print-directory -C Experiment/Forsterite pdgui-regenerate
	@$(MAKE) --no-print-directory -C Experiment/drude-lorentz pdgui-regenerate
	@$(MAKE) --no-print-directory -C Experiment/interpolation pdgui-regenerate
	@$(MAKE) --no-print-directory -C Experiment/Mayerhofer pdgui-regenerate
	@$(MAKE) --no-print-directory -C SingleCrystal/Bi2Se3 pdgui-regenerate
	@$(MAKE) --no-print-directory -C SingleCrystal/Bi2Se3_film pdgui-regenerate
