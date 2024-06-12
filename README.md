# Ecohydrology_MetaEcosysem

Code supporting the manuscript "Mechanistic insights on riverine meta-ecosystems: network shape drives spatial biodiversity and trophic structures" by Luca Carraro and Hsi-Cheng Ho, published in *Ecohydrology* (2024).

## Content

- `BUILD_FOODWEBS.R`: generates the different food-web realizations.
- `BUILD_OCN.R`: creates the elongated and compact OCNs and produces Figures 1, S1a and S7.
- `RUN_ME.m`: runs metaecosystem model for a given OCN and a given food-web realization. Results are stored in the `results` subfolder.
- `ANALYZE_DATA.m`: analyses output of the meta-ecosystem model and produces manuscript figures.
- `results`: contains output of the meta-ecosystem model (one file per OCN and food-web realization).
- `utilities`: contains utility functions and intermediate datasets (e.g., stored OCNs).
