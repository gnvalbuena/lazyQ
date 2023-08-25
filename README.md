# **lazyQ**

***An R package for easy calculation of qPCR results from QuantStudio output tables***

## Installation

The package can be installed in R by entering the following commands:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("gnvalbuena/lazyQ")
```

## Usage

1. Load package and import raw qPCR data

```r
library(lazyQ)
rawdata <- readQ("path/to/folder")
```

The `readQ()` function will take either an individual file path to a QuantStudio output spreadsheet, a list of file paths (as strings), or a folder containing QuantStudio output spreadsheets and import all spreadsheets. Check help with `?readQ()` for examples.

2. Check that samples and gene naming are correct

```r
importCheck(rawdata)
```

This will produce a table of samples and genes included in the imported raw data, so any errors in sample or gene naming can be identified and corrected if necessary.

3. Calculate dCT and ddCT for each target gene.

  3.1. Calculate dCT relative to measurements of the reference housekeeeping gene for each sample. The dCT values will be calculated for all target genes (i.e. genes other than those specified as housekeeping). This step is currently set up for use of only one reference housekeeping gene.

```r
dCT.table <- analyzeQ.CT(rawdata, housekeepinggene = "GAPDH")
```

  

  3.2. Calculate ddCT relative to mean dCT measurements of specified control samples.

```r
ddCT.table <- analyzeQ.ddCT(dCT.table, control.samples = c("Control_1", "Control_2", "Control_3"))
```

4. Alternatively, both steps in (3) can be performed as a single operation

```r
ddCT.table <- analyzeQ(readqtable, housekeepinggene = "GAPDH", control.samples = c("Control_1", "Control_2", "Control_3"))
```

