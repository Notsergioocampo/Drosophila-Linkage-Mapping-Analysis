# Drosophila Linkage Mapping Analysis

This repository contains a complete genetic linkage mapping analysis of three recessive mutations in *Drosophila melanogaster* (fruit flies). The project demonstrates classical three-point testcross analysis to determine gene order, calculate map distances, and construct genetic maps based on recombination frequencies.

## Project Overview

Genetic linkage mapping is a fundamental technique in genetics that uses recombination frequencies to determine the relative positions of genes on chromosomes. This project analyzes three mutant traits in fruit flies:

- **Crimson (cn)**: Affects eye color, producing deep red eyes
- **No Take Off (nto)**: Affects wing morphology, resulting in flightless wings  
- **Spray Tan (sp)**: Affects body pigmentation, creating darker body color

## Scientific Background

### Genetic Linkage Principles
When genes are located on the same chromosome, they tend to be inherited together rather than assorting independently. The frequency of recombination between linked genes is proportional to their physical distance, allowing researchers to construct genetic maps measured in centimorgans (cM).

### Three-Point Testcross Analysis
This powerful technique uses three linked genes to:
1. Determine gene order by analyzing double crossover events
2. Calculate precise map distances between loci
3. Detect genetic interference (suppression of nearby crossovers)
4. Compare experimental results with established genetic maps

## Project Structure

```
Drosophila-Linkage-Mapping-Analysis/
├── R/
│   └── genetic_mapping_functions.R    # Core analysis functions
├── reports/
│   ├── drosophila_report.Rmd         # Complete analysis report
│   ├── drosophila_report.html        # Rendered report output
│   └── references.bib                # Academic references
├── data/                             # Experimental data (referenced in report)
├── figs/                             # Generated figures and plots
├── Drosophila_Linkage_Mapping.Rproj  # RStudio project file
└── README.md                         # This file
```

## Key Findings

### Gene Order and Map Distances
- **Gene Order**: Crimson — No Take Off — Spray Tan
- **Crimson to No Take Off**: 32.5 cM
- **No Take Off to Spray Tan**: 26.8 cM
- **Total Map Length**: 59.4 cM

### Statistical Validation
- **Chi-Square Test**: χ² = 956.16, p < 0.0001
- **Conclusion**: Strong evidence against independent assortment, confirming genetic linkage
- **Interference**: I = 0.007 (minimal crossover suppression)

### Comparison with Standard Map
The calculated gene order and distances align closely with established *Drosophila* Chromosome 2 markers:
- Crimson corresponds to *brown* (bw) at ~104.5 cM
- No Take Off corresponds to *vestigial* (vg) at ~67.0 cM  
- Spray Tan corresponds to *black* (b) at ~48.5 cM

## Methodology

### Experimental Design
1. **Reciprocal Crosses**: Confirmed autosomal recessive inheritance
2. **F1 Testcross**: Generated recombinant progeny for mapping
3. **Phenotypic Scoring**: Classified 1,457 F2 offspring into 8 phenotypic classes
4. **Recombination Analysis**: Calculated crossover frequencies between loci

### Data Analysis Pipeline
1. **Phenotype Classification**: Categorized offspring as Parental, Single Crossover, or Double Crossover
2. **Recombination Frequency Calculation**: RF = (Recombinants/Total) × 100
3. **Map Distance Conversion**: 1% recombination = 1 centimorgan
4. **Statistical Testing**: Chi-square analysis to test for linkage
5. **Interference Analysis**: Coefficient of Coincidence calculation

## Core Functions

The `genetic_mapping_functions.R` file contains specialized functions for:

- **`build_triple_map()`**: Main analysis function that classifies phenotypes and calculates map distances
- **`draw_genetic_map()`**: Creates publication-quality genetic map visualizations
- **`calculate_chi_square()`**: Performs statistical testing for linkage
- **`calculate_interference()`**: Analyzes crossover interference patterns
- **`plot_phenotype_counts()`**: Visualizes phenotypic distributions
- **`calculate_map_confidence_intervals()`**: Provides statistical confidence for map distances

## Usage

To reproduce the analysis:

```r
# Load the analysis functions
source("R/genetic_mapping_functions.R")

# Input your phenotype data
df <- data.frame(
  phenotype = c("all wt", "eye-", "wing-", "body-", "eye- wing-", 
                "eye- body-", "body- wing-", "eye- wing- body-"),
  count = c(547, 233, 62, 148, 117, 64, 115, 171)
)

# Perform the complete analysis
map_results <- build_triple_map(df)

# View the genetic map
print(map_results$map_summary)

# Create visualizations
draw_genetic_map(map_results)
plot_phenotype_counts(map_results$data)
```

## Educational Value

This project demonstrates several fundamental genetic concepts:

1. **Mendelian Inheritance**: Autosomal recessive trait transmission
2. **Genetic Linkage**: Deviation from independent assortment
3. **Recombination**: Physical basis of genetic mapping
4. **Statistical Analysis**: Hypothesis testing in genetics
5. **Map Construction**: Converting recombination data to genetic distances
6. **Interference**: Crossover suppression mechanisms

## Applications

Linkage mapping has numerous applications in genetics and genomics:

- **Gene Discovery**: Identifying disease gene locations
- **Breeding Programs**: Marker-assisted selection in agriculture
- **Evolutionary Studies**: Comparative genome analysis
- **Medical Genetics**: Understanding inheritance patterns
- **Genomic Research**: Chromosome structure and organization

## Technical Requirements

- **R** (version 4.0 or higher)
- **RStudio** (recommended for report generation)
- Required R packages: `tidyverse`, `ggplot2`, `knitr`, `dplyr`

## Academic Context

This project was conducted as part of genetics coursework at the University of San Francisco. The analysis follows established protocols for *Drosophila* linkage mapping and provides students with hands-on experience in classical genetic analysis techniques.

## References

The complete analysis includes citations to foundational genetics literature, including works by Mendel, Morgan, Sturtevant, and modern genomic databases. See `reports/references.bib` for the complete bibliography.

## Citation

If you use this analysis or code in your research, please cite:

Ocampo, S. (2025). Drosophila Linkage Mapping Analysis. University of San Francisco, Department of Biology. GitHub repository: https://github.com/Notsergioocampo/Drosophila-Linkage-Mapping-Analysis

## License

This project is licensed under the MIT License - see the LICENSE file for details.
