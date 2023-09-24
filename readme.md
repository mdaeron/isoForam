[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8367191.svg)](https://doi.org/10.5281/zenodo.8367191)

# isoForam

Full code base and data for *Revisiting oxygen-18 and clumped isotopes in planktic and benthic foraminifera*, by Daëron & Gray (2023). Preprint available [here](http://daeron.fr/biblio/Daeron-Gray-2023.pdf).

Any questions/suggestions/comments/rants? Please open an issue on github or [contact the authors](mailto:daeron@lsce.ipsl.fr,william.gray@lsce.ipsl.fr?subject=[isoForam]).

## Install requirements

To install the required Python environemnt:

```
conda env create -f environment.yml
conda activate isoForam
```

It may be faster to use `mamba` instaed of `conda` in the first line:

```
mamba env create -f environment.yml
conda activate isoForam
```


## Executing the code

To execute the code, simply call `make` from the root directory. Most files output during the previous run will be deleted, so that each code run will start from a clean state. But this means that you should copy any output files you wish to save to another directory.

Throughout the code, some expensive calculations are skipped by default, using instead values stored after earlier runs. This should be obvious based on such boolean flags as:

```py
RUN_EXPENSIVE_CALCULATION = False
```

Simply change these flags to `True` to run the computation from scratch.

## Large data sets

The *Breitkreutz et al.* (2018), WOA23, and GLODAPv2 data are not included here but conveniently located `readme` files are included in their place, with the corresponding URLs and instructions as to which files are needed.

## CSV output

Most files generated in the code run will be copied to the `output` directory. One of those files, named `foram_D47_calibration_data.csv`, provides a very detailed listing of all properties computed for each sample in the clumped-isotope data set. Here is a list of the fields in this table:

* **`Sample`**: sample name, usually corresponding to a unique combination of `Site` & `Species`
* **`Species`**: foraminifer species
* **`Type`**: `benthic` or `planktic`
* **`Site`**: core-top site name
* **`Lat`**: core-top latitude
* **`Lon`**: core-top longitude
* **`Depth`**: for benthic samples, the core-top depth; for planktics, the assumed range of calcification depths
* **`d13C_VPDB`**: δ13C value of the sample in the VPDB scale
* **`SE_d13C_VPDB`**: SE of the above
* **`d18O_VPDB`**: δ18O value of the sample in the VPDB scale
* **`SE_d18O_VPDB`**: SE of the above
* **`D47_ICDES`**: Δ47 value of the sample in the I-CDES scale
* **`SE_D47_ICDES`**: SE of the above
* **`Tpub`**: for benthic samples, originally published calcification temperature
* **`SE_Tpub`**: SE of the above
* **`Twoa23`**: Calcification temperature estimate from atlas (WOA23); for planktics, based on assumed calcification depths
* **`SE_Twoa23`**: SE of the above
* **`Twoa23_500m`**: Calcification temperature estimate from atlas (WOA23); for planktics, based on the whole interval 0-500 m
* **`SE_Twoa23_500m`**: SE of the above
* **`Twoa23_1500m`**: Calcification temperature estimate from atlas (WOA23); for planktics, based on the whole interval 0-1500 m
* **`SE_Twoa23_1500m`**: SE of the above
* **`Tbottom_woa23`**: Bottom temperature estimate from atlas (WOA23); for benthics, identical to `Twoa23`
* **`SE_Tbottom_woa23`**: SE of the above
* **`Tk18`**: Calcification temperature estimate from Breitkreuz et al. (2018); for planktics, based on assumed calcification depths
* **`SE_Tk18`**: SE of the above
* **`Tk18_500m`**: Calcification temperature estimate from Breitkreuz et al. (2018); for planktics, based on the whole interval 0-500 m
* **`SE_Tk18_500m`**: SE of the above
* **`Tk18_1500m`**: Calcification temperature estimate from Breitkreuz et al. (2018); for planktics, based on the whole interval 0-1500 m
* **`SE_Tk18`**: SE of the above
* **`Tiso_species_offset`**: which B value (from table 4 in _Daëron & Gray_ in  was used to compute 18α
* **`Tiso_species`**: Calcification temperature estimate based on the 18α calibrations compiled in table 4 of _Daëron & Gray_; for planktics, based on assumed calcification depths
* **`SE_Tiso_species`**: SE of the above
* **`Tiso_species_500m`**: Calcification temperature estimate based on the 18α calibrations compiled in table 4 of _Daëron & Gray_; for planktics, based on the whole interval 0-500 m
* **`SE_Tiso_species_500m`**: SE of the above
* **`Tiso_species_1500m`**: Calcification temperature estimate based on the 18α calibrations compiled in table 4 of _Daëron & Gray_; for planktics, based on the whole interval 0-1500 m
* **`SE_Tiso_species_1500m`**: SE of the above
* **`Tiso_KON97`**: Calcification temperature estimate based on _Kim & O'Neil_ (1997); for planktics, based on assumed calcification depths
* **`SE_Tiso_KON97`**: SE of the above
* **`Twoa23_vs_Tiso_species`**: whether the sample is “discordant” or not, as defined by _Daëron & Gray_, based on assumed calcification depths
* **`Twoa23_vs_Tiso_species_500m`**: whether the sample is “discordant” or not (as defined by _Daëron & Gray_), based on the whole interval 0-500 m
* **`Twoa23_vs_Tiso_species_1500m`**: whether the sample is “discordant” or not (as defined by _Daëron & Gray_), based on the whole interval 0-1500 m
* **`Ref`**: Which studie this sample first appeared in
