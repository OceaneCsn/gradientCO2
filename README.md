# Context

We generated a transcriptomic dataset of the root response to gradually increasing CO2 concentrations in *Arabidopsis thaliana*. In addition to allowing a better resolution for GRN inference, measuring gene expression under a full range of CO2 concentrations has the potential to shed light on the dynamic of CO$_2$ response : is gene expression linearly reprogrammed as CO2 rises, or are there step functions and abrupt changes at specific CO$_2$ levels? Furthermore, we made the decision to investigate different types of N sources for the plant : nitrate and ammonium nitrate nutrition. This was motivated by the observation in the literature that ammonium nutrition elicits different phenotypic responses than nitrate nutrition in the face of CO$_2$ elevation. 


Arabidopsis Columbia ecotypes were hydroponically grown in 5 different controlled chambers, differing only in their CO$_2$ concentrations : 400, 525, 650, 775 and 900ppm. Inside a chamber, plants were separated in two groups, one receiving nitrate (KNO$_3$) and one receiving an ammonium nitrate (KNO$_3$-NH$_4$) mix, both resulting in a N concentration equal to 0.5 mM. Experiment carried out in the Ecotron in november 2021.


# Phénotypes analyses

Increase in biomass as CO$_2$ increases, even more marked in mix than KNO3. See `phenotypes/biomass.html` for graph and splines analyses.

Decrease in N content as CO$_2$ increases, even more marked in mix than KNO3. See `phenotypes/N-C.html` for graph and splines analyses.

Correlation analyses between N content and biomass are also done, and show no significant association for comparable N nutrition and CO2 condition.

# Transcriptomics analyses

A global differential expression analysis was carried out using splines modelling and edgeR, and the a GRN of the root response to a CO2 gradient was inferred using the integrative inference method : [bRF](https://github.com/OceaneCsn/integrative_GRN_N_induction). See `grn_inference/` folder for expression data preparation, TF binding sites data preparation, and inference scripts. The functions for bRF and its C++ dependencies and in `inference_functions/`.

# Results

Results in the form of figures for phentotypic analyses and GRN inference are stored in the folder `results/`.