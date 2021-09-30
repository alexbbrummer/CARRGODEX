# CARRGODEX
This repository contains the raw and processed data and code used for the analysis and figure generation of the manuscript Destabilization of CAR T-cell treatment efficacy in the presence of dexamethasone.

The files 2009011552P1R_Export.xlsx and 2009011553P1R_Export.xlsx are the raw outputs from the xCelligence cell killing assay experiments.  The first file corresponds to cell lines PBT030 (PBT1) and PBT128 (PBT2), while the second file corresponds to PBT138 with medium antigen levels (PBT3) and PBT138 with high antigen levels (PBT4).

The files car_T_dex_pbt1_pbt2_lines.csv and car_T_dex_pbt3_pbt4_lines.csv are processed data files consiting of the time-series values for average measurements and ranges of measurements for all treatment scenarios.

The file cytometry_results.xlsx containts the initial and final measurements of tumor and CAR T-cell populations for all treatment scenarios.

The Jupyter notebooks car_T_dex_fiting_pbt1.ipynb, car_T_dex_fiting_pbt2.ipynb, car_T_dex_fiting_pbt4.ipynb contain all of the algorithmic routines used to fit the modified CARRGO model to the experimental data.

The fit_parameters.csv contains all of the modified CARRGO model parameters, including standard errors and experimental conditions.

The Jupyter notebook graphing_fits.ipynb was used to generate the graphs that appear in Figures 3, 4, 6a, 6c, 6e, and Supplementary Information S1 Fig, S2 Fig, and S3 Fig.

The Jupyter notebook eigenvalue_plotting.ipynb was used to generate the graphs that appear in Figure 5.

The R file barplots_success_threshold_plotting.R was used to generate the barplots that appear in Figure 6b, 6d, 6f and the graph in Figure 7.

The file flow_ci_coreelation_data.csv contains measurements of xCelligence cell index and flow cytometry cell count used to establish correlation strength between these to indicators of tumor cell presence.  This data appears in Supplementary Information S1 File.

The R file flow_ci_correlation_graphing.R produces the graphical and tabular analysis of corrrelation strength between xCelligence cell index and flow cytometry cell count.  These results appear in Supplementary Information S1 File.
