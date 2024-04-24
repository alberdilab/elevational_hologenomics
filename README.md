![elevationhologenomics](https://github.com/alberdilab/elevational_hologenomics/assets/60347377/5a97f021-61d8-4757-9096-6815f61d8fe8)


# Elevational hologenomics

Data and code repository of the population-level metagenomic data of lizards across elevational gradients in various mountain ranges of the Pyrenees. The project aims to understand the role of gut microbes in the adaptation to different climatic regions using a hologenomic approach.

## 1. Study design

Four different mountain transects, two on the Spanish side and two on the French side, were used to sample lizards across altitudinal gradients ranging from 900m to 2400 ASL. Per sampling point 4-5 individuals were sampled to obtain faecal samples.


## 2. Bioinformatic procedures

Data processing to generate annotated metagenome-assembled genomes and genome count tables was conducted using the following Snakemake pipeline: [EHI_bioinformatics](https://github.com/earthhologenome/EHI_bioinformatics/tree/mjolnir). Data analysis procedures source from the outputs of this pipeline.


## 3. Data analysis

Samples from the four transects were analysed together and the code used can be found in the Rmd files stored in the root directory of this repository, while the bookdown-rendered webbook is available at:

[alberdilab.github.io/elevational_hologenomics](https://alberdilab.github.io/elevational_hologenomics/)

While the webbook provides a user-friendly overview of the procedures, analyses can be directly reproduced using the Rmd documents. Note that the code chunks that require heavy computation have been tuned off using 'eval=FALSE'. To re-render the webbook, you can use the following code:

```
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")
```
