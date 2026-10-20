# corto 1.4.0

## New features

* `corto()` gains a `dpi` argument to select the Data Processing Inequality
  rule. The default, `"target"`, is the rule used so far: for each target only
  the edge from the centroid with the strongest absolute correlation survives,
  and results are identical to previous versions. `"triangle"` applies the
  ARACNe rule, removing an edge only if it is the weakest in a
  centroid-centroid-target triangle; it produces larger networks, closer to
  those of ARACNe, at a considerably higher computing time. Proposed and
  implemented by Hugo Tovar (#16).

* New exported `dpi_edges()`, which applies either DPI rule to any list of
  centroid-target edges.

## Documentation

* The package description, the vignette and the `corto()` help page now
  describe the DPI step as it is actually performed. They previously
  described the ARACNe triangle rule, while `corto()` applied the per-target
  rule (#16).
* The vignette no longer states that correlations are calculated on ranks:
  `corto()` uses Pearson correlation on the input matrix as provided. Rank the
  matrix beforehand to work on Spearman correlations.

# corto 1.3.1

## New features

* New `filter_regulon()` and `getregulon()` helpers for post-processing a
  regulon: subset it by centroid, likelihood and absolute correlation, and
  flatten it into a data frame of edges with optional tab-separated export.
  Proposed by Hugo Tovar (#15).

* `mra()` now performs Signature Master Regulator Analysis when `expmat1` is
  provided as a named vector, as the documentation had always described. The
  null model permutes the signature values across feature names, and `nperm`
  defaults to 1000 in this mode. Thanks to Hugo Tovar for the detailed
  diagnosis and a reference implementation (#13).

## Breaking changes

* Sample-by-sample `mra()` (a single `expmat1`, no `expmat2`) now returns the
  same list as the other modes, `list(nes, pvalue, sig, regulon)`, instead of a
  bare matrix. Use `mra(expmat, regulon=regulon)$nes` to get the previous
  output (#1).

## Bug fixes

* `corto()` failed with `invalid 'row.names' length` when a single centroid was
  provided. The internal correlation step dropped to a vector and lost its row
  names. Thanks to Hualin Wang for the fix (#6, #14).
* Sample-by-sample `mra()` returned an all-NA matrix whenever any regulon target
  had zero variance in the input matrix. The permuted null signatures were not
  NA-guarded the way the real signature was.
* `gsea()` errored on R >= 4.2 when `method` was left at its default, since the
  default is a length-2 vector. It now uses `match.arg()` and defaults to
  `"permutation"` as before.
* `plot_gsea(omit_middle=TRUE)` errored on an undefined `legend_position`.
* `mraplot()` errored on regulons with fewer than 12 targets, which is reachable
  with the default `minsize=10`.
* `ssgsea(scale=TRUE)` silently dropped the sample names from the returned NES
  matrix.
* `mra()` failed on regulons left with a single centroid after `minsize`
  filtering.
* `mraplot()` now stops with an informative message when given sample-by-sample
  results, which it cannot plot.
* The CNV correction step dropped to a vector when `cnvmat` and `inmat` had a
  single sample or a single target in common, the same dimension-drop problem
  fixed elsewhere. Guarded with `drop=FALSE`.
* `corto()` now stops with an informative message when no edge passes the
  correlation threshold, instead of failing on an invalid row name.

## Performance

* `corto()` is roughly 1.2x faster single-threaded and 1.4x faster on 4 threads.
  The input matrix is transposed once instead of once per bootstrap and is no
  longer shipped twice to the workers, DPI selection is done by a radix sort
  rather than a grouped data frame, and the regulon is assembled in one pass
  instead of rescanning the edge table for every centroid. Results are
  bit-identical to previous versions.

## Other

* `dplyr` and `gplots` are no longer required, which removes about 20 recursive
  dependencies and makes installation considerably faster. `knitr` and
  `rmarkdown` moved from Imports to Suggests, `grDevices` and `graphics` added.
* Documented defaults corrected for `scatter(bgcol=)`, `val2col(nbreaks=)` and
  `plot_gsea2()`.
