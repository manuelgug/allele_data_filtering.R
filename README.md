# allele_data_filtering.R

## Description

`allele_data_filtering.R` applies contaminants filtering and minimum allele frequency (MAF) filtering to remove potential false positives from [mad4hatter's](https://github.com/EPPIcenter/mad4hatter) main output allele_data.txt and resmarker_microhap_table.txt.

## Dependencies

- `dplyr`
- `tidyr`
- `optparse`
- `ggplot2`
- `gridExtra`

## Usage

```bash
Rscript allele_data_filtering.R [--allele_table PATH] [--resmarkers_table PATH] [--CFilteringMethod METHOD] [--MAF VALUE] [--exclude_file PATH]
```

- `--allele_table`: Path to the input allele table.

- `--resmarkers_table`: Path to the input resmarkers table (optional).

- `--CFilteringMethod`: Contaminants filtering method. Options: "global_max", "global_q95", "amp_max", "amp_q95". Default: "global_max".

  - "global_max": Single threshold derived from the maximum read count from all amplicons across negative controls.
  - "global_q95": Single threshold derived from the 95th percentile of read counts from all amplicons across negative controls.
  - "amp_max": Amplicon-specific thresholds derived from the maximum read count for each amplicon across negative controls.
  - "amp_q95": Amplicon-specific thresholds derived from the 95th percentile of read counts for each amplicon across negative controls.

- `--MAF`: Minimum allele frequency filter. Default: 0.

- `--exclude_file`: Path to the file containing sampleIDs to exclude (optional).

## Example

```bash
Rscript allele_data_filtering.R --allele_table allele_data.txt --resmarkers_table resmarker_microhap_table.txt --CFilteringMethod amp_max --MAF 0.01 --exclude_file samples_to_exclude.txt
```
### Report visualization
![report_visualization.jpg](https://github.com/manuelgug/allele_data_filtering.R/blob/main/report_visualization.png)

## Nomenclature of controls

- **Positive Controls**: The script identifies positive controls using the following conditions:
  - Sample IDs containing "3D7" (case-insensitive) and not containing "Dd2", "HB3", or "PM" (case-insensitive).

- **Negative Controls**: The script identifies negative controls using the following conditions:
  - Sample IDs containing "Neg" (case-insensitive).
