# Figure reproduction

`00_Run_all_figures.R` contains the integrated workflow for main Figures 1–4 and Supplementary Figures 1–4. It starts with an annotated Seurat object, not raw data. Upstream preprocessing, QC, normalization and annotation scripts are not included in this draft.

Install the required R packages: Seurat, data.table, ggplot2, patchwork, scales, FNN and Matrix. For image reading, install magick or the appropriate png, jpeg or tiff package. Actual package versions must be recorded from the original analysis environment.

Place the annotated object at `6 Analysis/Data objects/CRC_CosMx_final_annotated_full_object.rds` within your local project directory. Set `CRC_COSMX_SAMPLE` to the sample identifier in the approved object. The public default is a placeholder.

```r
Sys.setenv(CRC_COSMX_PROJECT_DIR = "/path/to/local/project", CRC_COSMX_SAMPLE = "PUBLIC_SAMPLE_ID")
source("scripts/00_Run_all_figures.R")
```

For Figure 1, supply the approved CellComposite images in `scripts/Figure1_assets/` as `Figure1_CellComposite_Non-response.jpg` and `Figure1_CellComposite_Response.jpg`, or configure the image paths in the script. Missing images may result in placeholders. Figures and tables are written beneath `6 Analysis/Figures` and `6 Analysis/Tables` of the local project. The script also saves session information.

This draft has not been executed against the original dataset and requires author verification before release.
