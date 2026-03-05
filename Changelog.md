**v1.3 Changes**

- PCA is calculated using a dynamic number of HVGs (specified by asymptote detection) as well as using ALL genes
- Switched plotting over to plotnine (updated pcaplot_env.yaml)
- Added fast clustering (update pcaplot_env.yaml)
- Added log file output that informs how many HVGs were used for PCA calculation
- Updates to multiqc_env.yaml to help github tests pass
