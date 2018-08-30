# chikv-colombia-emod
Code for model of chikungunya in Colombia by Moore et al. 2018 published in BMC Medicine https://doi.org/10.1186/s12916-018-1127-2

Repository content:

Eradication folder contains the EMOD v2.10 executable with the DTK-Dengue extension. See https://github.com/InstituteforDiseaseModeling/EMOD for model documentation and latest version of the executable.

Eradication/InputDataFiles folder contains climate, demographics, and human movement files for each department in Colombia.

The dtk-tools-fork folder contains the package of python scripts used to manage EMOD simulations and calibrations. For details on the structure of DTK and the latest version of dtk-tools see https://github.com/InstituteforDiseaseModeling

The dtk-submissions folder contains python scripts used to run the various simulations and calibrations described in Moore et al. 2018.

The inputs folder contains the department-level time series of weekly reported chikungunya cases from the start of the epidemic through week 2 of 2016.

The R folder contains R scripts used to analyze model output.
