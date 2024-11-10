# HighEnergyGammarayDetectors

This is the code pipeline for processing the spectra obtained in the space detector lab. 
The SpaceDL_code file plots the spectra of the sources, fits a gaussian profile to a spectrum and prints the parameters of the profile. The centroid channel, fwhm, absolute efficiency and all associate errors are all obtained and added to a seperate table, which can be seen in the appendix of the report. The code that was changing was the detector, source and the mask, for the region of interest.
The SpaceDL_plot file plots the channel-, resolution- and efficiency-energy plots for the values obtained from the previous file. 
