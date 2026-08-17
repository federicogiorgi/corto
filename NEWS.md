# corto 1.2.2 (Development version)

## New features

- Added `filter_regulon()` to subset a regulon object by transcription factors, correlation strength, and likelihood threshold.
- Added `getregulon()` to convert a regulon object into a tidy tibble, optionally exporting to TSV.
- Included an example regulon in `inst/extdata/example_regulon.rds`.
- Added vignette `regulon_utils.Rmd` demonstrating how to filter and extract a regulon.
- Added runnable example script `inst/examples/regulon_utils_example.R`.

These utilities improve downstream analysis and visualization of regulons inferred using `corto()`.

