all: cleanup processD47 processd18O assignd18Osw Tcalcif Tatlas Tcalcif_onemoretime swchem finalplots malevich meckler d13C finalcleanup
	@echo "\nALL DONE."

cleanup:

	@echo "\nCLEANING UP...\n"

	rm -f 1_compile_D47_data/foram_D47_compilation.csv
	rm -f 1_compile_D47_data/*.pdf

	rm -f 2_compile_d18O_data/*.csv
	rm -f 2_compile_d18O_data/*.pdf
	rm -f 2_compile_d18O_data/plots/*/*.pdf

	rm -f 5_assign_d18Osw/*.csv
	rm -f 5_assign_d18Osw/plots/*/*.pdf

	rm -f 6_assign_atlas_T/*.csv
	rm -f 6_assign_atlas_T/histograms/pdf/*.pdf

	rm -f 7_assign_Tcalcif/*.csv
	rm -f 7_assign_Tcalcif/*.pdf

	rm -f 8_seawater_chemistry/*.csv

	rm -f 9_final_plots/sites.csv
	rm -f 9_final_plots/plots/*.pdf

	rm -f output/*.py
	rm -f output/*.pdf
	rm -f output/*.csv

processD47:

	@echo "\nPROCESSING D47 DATA..."
	@echo; cd 1_compile_D47_data; ./plot_as_published.py; ./foram_D47_compilation.py
	
processd18O:
	@echo "\nPROCESSING d18O DATA..."
	@echo; cd 2_compile_d18O_data; ./foram_d18O_compilation.py; ./species_offsets.py

assignd18Osw:
	@echo "\nEXTRACTING d18Osw DATA FROM ATLAS..."
	@echo; cd 5_assign_d18Osw; ./estimate_bottom_d18Osw.py; ./estimate_planktic_d18Osw.py

Tcalcif:
	@echo "\nCOMPUTE CALCIFICATION TEMPERATURES..."
	@echo; cd 7_assign_Tcalcif; ./assign_Tcalcif.py

Tatlas:
	@echo "\nEXTRACTING TEMPERATURE DATA FROM ATLAS..."
	@echo; cd 6_assign_atlas_T; ./assign_atlas_T.py

Tcalcif_onemoretime:
	@echo "\nCOMPUTE CALCIFICATION TEMPERATURES AGAIN..."
	@echo; cd 7_assign_Tcalcif; ./assign_Tcalcif.py

swchem:
	@echo "\nSEAWATER CHEMISTRY..."
	@echo; cd 8_seawater_chemistry; ./seawater_chemistry.py

finalplots:
	@echo "\nCREATE FINAL PLOTS..."
	@echo; cd 9_final_plots; ./Tiso_vs_Tiso_500m.py; ./compile_sites.py; ./sitemap.py; ./regression_plots.py

malevich:
	@echo "\nMALEVICH (2019)..."
	@cd Malevich_2019; ./4_plots.py

meckler:
	@echo "\nMECKLER (2022)..."
	@echo; cd Meckler_2022; ./meckler_2022.py

d13C:
	@echo "\nCHECK d13C..."
	@echo; cd check_d13C; ./check_d13C.py

finalcleanup:
	@echo "\nCLEAN UP..."
	@cp 1_compile_D47_data/*.pdf ./output/
	@cp 2_compile_d18O_data/*.pdf output/
	@cp 2_compile_d18O_data/foram_d18O_compilation.csv output/
	@cp 2_compile_d18O_data/species_offsets.csv output/
	@cp 5_assign_d18Osw/*.pdf output/
	@cp 7_assign_Tcalcif/*.pdf output/
	@cp 7_assign_Tcalcif/foram_D47_calibration_data.csv output/
	@cp 9_final_plots/plots/*.pdf output/
	@cp 9_final_plots/sites.csv output/
