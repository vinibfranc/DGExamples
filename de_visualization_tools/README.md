# Visualization tools tests

## Tools installed inside R:

1) (Failing)[PIVOT](https://github.com/qinzhu/PIVOT)
    * Command line dependencies: <br>
        ```sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libcurl4-gnutls-dev libcurl4-openssl-dev libcairo2-dev```
    * R packages: <br>
    ```
    # dependecies that needs to be manually installed 
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install()
    BiocManager::install("RCurl")
    install.packages("devtools") 

    library("devtools")
    source("http://bioconductor.org/biocLite.R")  
    biocLite("GO.db")
    biocLite("HSMMSingleCell")
    biocLite("org.Mm.eg.db")
    biocLite("org.Hs.eg.db")
    biocLite("DESeq2")
    biocLite("SingleCellExperiment")
    biocLite("scater")

    # Install PIVOT
    install_github("qinzhu/PIVOT")
    biocLite("BiocGenerics") # You need the latest BiocGenerics >=0.23.3
    ```

2) [Glimma](https://github.com/Shians/Glimma)
    * R packages:
    ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("Glimma")
    ```

3) [DEBrowser](https://github.com/UMMS-Biocore/debrowser)
    * R packages:
    ```
    if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
    BiocManager::install("debrowser")
    ```

4) [ExpressionDB](https://github.com/5c077/ExpressionDB)
    * Download and unzip file from Github repo and run global.R Shinny app!