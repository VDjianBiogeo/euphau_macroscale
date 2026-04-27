# Biogeographic partitioning of euphausiids across Southern Ocean fronts revealed by environmental partitioning models

*Purpose*\
This project aims to test the hypothesis that ACC fronts partition euphausiid niches through distinct environmental mechanisms using ensemble species distribution models. We quantify the relative contribution of environmental structures to the spatial organisation of 29 euphausiid species across the Southern Indian Ocean and the Indian sector of the Southern Ocean.

## Content

The repository of this project is structured as follow:

```         

.
├─ README.md                          # Description of the project

|

├─ analyses                           # All R scripts used for analysis
| ├─ mpi scripts                      # All R and shell scripts used in the SACADO Cluster for modelling
| ├─ raw_euphau_data_treatment.R      # R script for raw data extraction and treatment
| ├─ exploratory_analyses.R           # R script to explore the treated data's distribution and quality
| ├─ pre_mod_analyses.R               # R script to select environmental predictors and optimise algorithm hyperparameters
| ├─ post_mod_analyses.R              # R script to treat and analyse modelling outputs
|
|
|
├─ R/                                 # R function directory
|
|
├─ make.R                             # R script to launch the complete analysis pipeline
```

## Data

*Biological data*\
Euphausiid species observations were extracted at a worldwide scale from the Ocean Biodiversity Information System ([OBIS](https://obis.org)) and the Global Biodiversity Information Facility ([GBIF](https://www.gbif.org/fr/)) databases. To this pool of occurrences were added data from regional surveys, such as the REPCCOAI program (Djian et al., 2025; Merland et al., 2025). Data were extracted at a worldwide scale between 1950 to time of download (13/11/2023). Raw data are available on [Zenodo](10.5281/zenodo.19815422) in the input/raw/bio directory.

*Environmental predictors*\
“Physical” variables were extracted from the [Global Ocean Physics Reanalysis](https://doi.org/10.48670/moi-00021) available on [Copernicus](https://marine.copernicus.eu/?pk_vid=b2a6239971064277177185355142c261) in october of 2023. “Biogeochemical” variables were extracted from the [Global Ocean Biogeochemistry Hindcast](https://doi.org/10.48670/moi-00019) on [Copernicus](https://marine.copernicus.eu/?pk_vid=b2a6239971064277177185355142c261) in october of 2023. Raw data are available on [Zenodo](10.5281/zenodo.19815422) in the input/raw/env directory.

## Workflow

After spatially joining every euphausiid species occurrence to environmental predictors, we ran a two-step environmental predictors pre-selection procedure in order to limit multicollinearity. Three pseudo-absences datasets were generated using surface range enveloppe. Hyperparameter optimisation for Random Forest and Gradient Boosted Models was done using the 'h2o' R package (Fryda et al., 2024) for each pseudo-absence dataset. Models were calibrated using a five-fold cross validation repeated 5 times for each pseudo-absence dataset, resulting in 150 simulations per species. Each simulation was evaluated using the AUC metric. All simulations with an AUC higher than 0.7 were then stacked into an ensemble model.

![Overview of the modelling workflow implemented in this study, including pseudo-absences (PA) generation, calibration of individual models, construction of ensemble models, model evaluation, and spatial projection. AUC: Area Under the receiver operating characteristic Curve.](./Figure%20methodological%20framework%20axe%201%20v4.jpg)

In addition to full ensemble models including all selected predictors, three partial models were constructed for each species using only predictors associated with one environmental category: *Water mass properties*, *Oceanic dynamics*, or *Biological production*. Partial models were calibrated using the same modelling procedure as full models, differing only in the predictor set. Spatial projections were compared using semi-partial Spearman's rank correlation coefficients computed across the study area. A Principal Component Analysis was applied on the resulting matrix of semi-partial correlation coefficients to position species according to the environmental structure most closely approximating the full model. A hierarchical ascendant classification was performed on species coordinates in the PCA factorial space to identify species sharing similar environmental dependency profiles.

![Schematic representation of the environmental partitioning framework. Complete and partial models were compared using semi-partial correlations to quantify the relative contribution of environmental variable groups, followed by multivariate analyses.](./Figure%20environmental%20partitioning.jpg)

All the analysis pipelines can be launched from the `make.R` file. You need to launch:

1) `source(here::here("analyses", "raw_euphau_data_treatment.R"))`;

2) `exploratory_analyses.R` pipeline;

3) `pre_mod_analysis` pipeline down to the `mod_prep` function;

4) Launch the scripts in the mpi scripts directory up to `6_ensemble_proj.R` (through cluster if needed);

5) Launch the `pre_mod_analysis` pipeline once the outputs from step 4 are available;

6) Launch the rest of the scripts in the mpi directory;

7) `post_mod_analysis.R`

## License

This project is released under the [GPL-2](https://choosealicense.com/licenses/gpl-2.0/) license.

## References

Djian, V., Merland, C., Thellier, M., Leroy, B., Cotte, C., & Koubbi, P. (2025). Species assemblages of key prey taxa from the Southern Ocean and the Southern Indian Ocean using biogeographic networks. *CCAMLR SCIENCE*, *25*, 19‑38.

Fryda, T., LeDell, E., Gill, N., Aiello, S., Fu, A., Candel, A., Click, C., Kraljevic, T., Nykodym, T., Aboyoun, P., Kurka, M., Malohlava, M., Poirier, S., Wong, W., Rehak, L., Eckstrand, E., Hill, B., Vidrio, S., Jadhawani, S., … H2O.ai. (2024). *h2o : R Interface for the « H2O » Scalable Machine Learning Platform* (Version 3.44.0.3) [Logiciel]. <https://cran.r-project.org/web/packages/h2o/index.html>

Merland, C., Thellier, M., Djian, V., Cotte, C., & Koubbi, P. (2025). Macrozooplankton from Crozet to Kerguelen and Subtropical Southern Indian Ocean. *CCAMLR Science*, *25*, 99‑116.
