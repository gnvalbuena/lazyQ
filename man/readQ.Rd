% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readQ.R
\name{readQ}
\alias{readQ}
\title{Import QuantStudio output files}
\usage{
readQ(x = "/path/to/QuantStudio.xls")
readQ(x = list("/path/to/QuantStudio1.xls", "/path/to/QuantStudio2.xls", "/path/to/QuantStudio3.xls"))
readQ(x = "/path/to/folder/")
}
\arguments{
\item{x}{Any of the following options: (1) a path to a single QuantStudio output file, (2) a list of paths to a QuantStudio output file, or (3) a path to a folder containing QuantStudio output files to be analyzed.}
}
\value{
A table of raw qPCR data from all imported QuantStudio output files
}
\description{
This function imports raw qPCR data from a single QuantStudio output file path,
a list of QuantStudio output file paths, or a path to a folder containing QuantStudio output files
}
