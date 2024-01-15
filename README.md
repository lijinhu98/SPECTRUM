# SPECTRUM

## Installation Guide
**Installing SPECTRUM**  
To install SPECTRUM, run:
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("lijinhu98/SPECTRUM")
```

## Tutorial
### Library
```
library(tidyverse)
library(reshape2)
```
### predict cell type
For query set matrix, rows should be cells and column should be genes (TPM).
```
cell.type <- celltype_predict(matrix)
cell.type
```
### predict patient group from sub cell type proportion
For query set proportion, rows should be sampleID and column should be proportion of each sub cell type in each person.
```
group.predict <- patient_group_predict(proportion)
group.predict
```
### predict patient group from cell X gene matrix
For query set matrix, rows should be cells and column should be genes (TPM) but the last column should be "sampleID" for each cell
```
sampleID <- matrix$sampleID
matrix <- matrix[,which(colnames(matrix) != "sampleID")]
group.predict <- merge_predict(matrix,sampleID)
group.predict
```

## Example
### Library
```
library(tidyverse)
library(reshape2)
```
### load data
```
load(file = "data/Qian.RData")
sampleID <- Qian$sampleID
Qian <- Qian[,which(colnames(Qian) != "sampleID")]
```
### predict cell type
```
cell.type <- celltype_predict(Qian)
cell.type
```
### predict patient group from cell X gene matrix
```
group.predict <- merge_predict(Qian,sampleID)
group.predict
```
