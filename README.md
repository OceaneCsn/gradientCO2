# Context


Experiment to quantify the effect of increasing CO$_2$ levels on *Arabidopsis thaliana*'s genomics and phenotypes, under two type of N nutrition. Carried out in the Ecotron in november 2021 by A. and O..

For 5 CO$_2$ levels, (400, 525, 650, 775, 900 ppm) and 2 N nutrition (KNO3 and mix NH4KNO3, both resulting in 0.5mM N availability), roots were collected and flash frozen, and leaves were partly flash frozen and partly dried for biomass and mineral status analyses.

Samples of roots ans shoots were samples by A. and L., milled by A. and O., RNAs extracted by A., aliquoted by A. 
RNA-Seq was performed by Genoscreen.

Samples of shoots to be dried were collected by O. and J.

# Phénotypes analyses

Increase in biomass as CO$_2$ increases, even more marked in mix than KNO3. See `phenotypes/biomass.html` for graph and splines analyses.

Decrease in N content as CO$_2$ increases, even more marked in mix than KNO3. See `phenotypes/N-C.html` for graph and splines analyses.

# Transcriptomics analyses

## Bioinformatics

Samples were processed on the plasticity server using the script `bioinformatics/batch_process.py` (open it for details and arguments used for the quality and mapping jobs). Quantification was done using `bioinformatics/quantif.sh` (open it to see the arguments passed to htseq-count).

High quality samples, high quality alignment rate, high fraction of reads quantified.


## Sample mix-ups are highly highly likely

Using PCA, corplot (`bioinformatics/format_raw_expression.R`) and NRT2.1 expression, I saw that replicates were not homogenous and that some N and CO2 conditions had been switched inadvertantly. I manually re-annotated the samples to have clean PCA, corplot showing maximum correlation between replicates, and coherent NRT2.1 expression profiles (high expression in KNO3, low in Mix). After discussion with A., we chose to use my corrected annotation as this is the most likely and coherent annotation (even though it is upsetting that we'll never be sure of the true annotation).

> Our RNA extraction must be clean, as A. made sucessful rt-QPCR with coherent profiles, so the mix-up must have been with the aliquots, or during the processing by Genoscreen...


## Statistical analyses

Splines analyses are envisioned to study the joint effect of CO2 and N nutrition on differential expression.
Moanin could also be used as the design is much more suited.