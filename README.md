# allele_data_filtering.R

## Description

`allele_data_filtering.R` applies contaminants filtering and minimum allele frequency (MAF) filtering to remove potential false positives from [mad4hatter's](https://github.com/EPPIcenter/mad4hatter) main output allele_data.txt and resmarkers_microhap_table.txt.

## Dependencies

- `dplyr` package

## Usage

```bash
Rscript allele_data_filtering.R <allele_table> <resmarkers_table> <CFilteringMethod> <MAF> <exclude_file>
```

- `allele_table`: Path to the input allele file.
- `resmarkers_table`: Path to the input resmarkers file.
- `CFilteringMethod`: Filtering method for contaminants. Options: `global_max`, `global_q95`, `amp_max`, `amp_q95`.

  - *global_max*: single threshold derived from the maximum read count from all amplicons across negative controls

  - *global_q95*: single threshold derived from the 95th percentile of read counts from all amplicons across negative controls

  - *amp_max*: amplicon-specific thresholds derived from maximum read count for each amplicon across negative controls

  - *amp_q95*:  amplicon-specific thresholds derived from 95th percentile of read counts for each amplicon across negative controls
  
- `MAF`: Minimum allele frequency. Default: 0.02.

- `exclude_file`: file with a list of samples to exclude (optional).

## Example

```bash
Rscript script_name.R allele_data.txt resmarkers_microhap_table.txt max 0.01 neg_controls_to_exclude.txt
```

## Nomenclature of controls

- **Positive Controls**: The script identifies positive controls using the following conditions:
  - Sample IDs containing "3D7" (case-insensitive) and not containing "Dd2", "HB3", or "PM" (case-insensitive).

- **Negative Controls**: The script identifies negative controls using the following conditions:
  - Sample IDs containing "Neg" (case-insensitive).
