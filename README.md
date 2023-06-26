# allele_data_filtering.R

## Description

`allele_data_filtering.R` applies contaminants filtering and minimum allele frequency (MAF) filtering to remove potential false positives from [mad4hatter's](https://github.com/EPPIcenter/mad4hatter) main output allele_data.txt.The filtered allele data is then exported to a file named "allele_data_filtered.txt".

## Dependencies

- `dplyr` package

## Usage

```bash
Rscript allele_data_filtering.R <path_to_file> <CFilteringMethod> <MAF>
```

- `path_to_file`: Path to the input allele data file.
- `CFilteringMethod`: Filtering method for contaminants. Options: `max`, `q95`, `amp_max`, `amp_q95`.

  - *max*: single threshold derived from the maximum read count from all amplicons across negative controls

  - *q95*: single threshold derived from the 95th percentile of read counts from all amplicons across negative controls

  - *amp_max*: amplicon-specific thresholds derived from maximum read count for each amplicon across negative controls

  - *amp_q95*:  amplicon-specific thresholds derived from 95th percentile of read counts for each amplicon across negative controls
  
- `MAF`: Minimum allele frequency. Default: 0.02.

## Example

```bash
Rscript script_name.R allele_data.txt max 0.01
```

## Control Nomenclature

- **Positive Controls**: The script identifies positive controls using the following conditions:
  - Sample IDs containing "3D7" (case-insensitive) and not containing "Dd2", "HB3", or "PM" (case-insensitive).

- **Negative Controls**: The script identifies negative controls using the following conditions:
  - Sample IDs containing "Neg" (case-insensitive).
