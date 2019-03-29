# 2019 cell line nicastrin knockdown paper
This repository contains the code required to generate figures and perform analysis for our paper knocking down nicastrin in two human cell lines.

## Related links
[GEO expression data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57949)
[Sequencing Project](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA268374)
[Miscellaneous data](https://figshare.com/projects/2019_Roberson_lab_nicastrin_cell_line_knockdown_data/58658)
[Preprint](https://www.biorxiv.org/content/early/2018/06/22/353482.1)
[Published paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/ced.13906)

## Requirements

### System
* R
* R Studio

### R libraries
* beadarray
* broom
* cowplot
* dplyr
* ggplot2
* ggrepel
* here
* limma
* pheatmap
* plyr
* purrr
* RColorBrewer
* readr
* reshape2
* stringr
* tibble
* tidyr
* UpSetR

## Running analysis

### Download data and run differential expression
1. Install R and R Studio if needed.
2. Clone this repository.
3. Install required R packages. Most of the tidyverse packages are CRAN hosted. Some others (beadarray, limma, etc) are Bioconductor packages.
4. Knit R Markdown files 00 - 03 in order (in the code directory). The 00 file will download the necessary files from FigShare. It will also auto-generate required directories. Output should save to the auto-generated "results" directory.

### Run Gene Ontology analysis
5. Using your browswer, navigate to the [gProfileR web interface](https://biit.cs.ut.ee/gprofiler).
6. Make sure the following options are checked:
 * Organism - Homo sapiens
 * Significant query
 * Ordered query
 * No electronic GO annotations
 * Hierarchical sorting - best parent moderate
 * Output type - textual (download)
 * User p-value - 0.05
 * Size of functional category Max - 500
 * Significance threshold - g:SCS
7. Profile the DE genes from HEK001 / HEK293 and increased / decreased separately by pasting them into the Query box.
8. Each profile will yield a text result. The filename will be something like "gprofiler_results_xxxxxxxxxxxxx.txt", where the X's are numbers. Save the file by adding the cell line and direction to the name, i.e. "gprofiler_results_xxxxxxxxxxxxx_hek001_down_v00.txt" or "gprofiler_results_xxxxxxxxxxxxx_hek001_up_v00.txt". Save the files in the auto-generated "data/gprofiler" directory in the repository. These files are required to generate the Gene Ontology enrichment figures. **Note** the original gprofiler output is included in the git repository.

### Generate remainder of analysis
9. Knit R Markdown files 04 - 15.
