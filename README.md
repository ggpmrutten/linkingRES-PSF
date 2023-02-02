# Analysis MS review Plant&Soil
by Gemma Rutten (gemma.rutten@unibe.ch)

Last edited 01/01/2023
* [preprint](https://europepmc.org/article/PPR/PPR517673) DOI: 10.21203/rs.3.rs-1813185/v1

### Investigates if published data on plant soil feedbacks can be explained by belowground root strategies

1. SpeciesPSF.R script collects data from publications (using A. Crawford.csv, B. Bennett.csv, C. Petermann.csv)
 * |- Prepares species lists for root trait database (sRoot_extraction.RES.PSF.R) see B below
 * |- Consistent formating of the plant soil feedback (PSF) variables (Linking PSF.RESpetben.R)
 
 * 1b. sRoot_extraction.RES.PSF.R collects root traits from database (sRoot.trait.Rdata) based on [GRoot](https://groot-database.github.io/GRooT/) 
  for species lists (study B and C) from script 1
 * |- Collection of root traits for the species lists
 * |- Calculation of the root economics spectrum

2. Linking PSF.RESpetben.R   
 Combines root traits (output sRoot_extraction.RES.PSF.R) with PSF measures (1.SpeciesPSF.R) 
 and calculates pairwise differences of PSFs in Root Economics Spectrum 

3. Tests.Figures.R stats and graphs based on output script 2 (PSF.RES.benpet.Rdata)

