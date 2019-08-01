# Visualization tools tests

## Tools installed inside R:

1) [PIVOT](https://github.com/qinzhu/PIVOT)
    * Command line dependencies: <br>
        ```sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libcurl4-gnutls-dev libcurl4-openssl-dev libcairo2-dev libxt-dev```
    * R packages: <br>
    ```
    # dependecies that needs to be manually installed 
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install()
    BiocManager::install("RCurl")
    install.packages("devtools") 

    library("devtools")
     
    BiocManager::install("GO.db")
    BiocManager::install("HSMMSingleCell")
    BiocManager::install("org.Mm.eg.db")
    BiocManager::install("org.Hs.eg.db")
    BiocManager::install("DESeq2")
    BiocManager::install("SingleCellExperiment")
    BiocManager::install("scater")
    BiocManager::install("BiocGenerics") # You need the latest BiocGenerics >=0.23.3

    # Install PIVOT
    devtools::install_github("qinzhu/PIVOT")
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