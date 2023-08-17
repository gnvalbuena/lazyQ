#' Import a single QuantStudio output file
#'
#' This function imports raw qPCR data from a QuantStudio output file
#'
#' @usage readQ.file(x = "/path/to/QuantStudio.xls")
#' @param x A path to a QuantStudio output file (i.e. `"/path/to/Quantstudio.xls"`)
#' @return A table of raw qPCR data from a QuantStudio output file
#'
#' @export
readQ.file <- function(x){

  filebase <- basename(x)

  outtable <- data.frame(file = filebase,
                         readxl::read_excel(path = x, sheet = "Results", skip = 46))

  outtable$CT[outtable$CT == "Undetermined"] <- NA
  outtable$CT <- as.numeric(outtable$CT)

  return(outtable)

}

#' Import QuantStudio output files from a list of files
#'
#' This function imports raw qPCR data from QuantStudio output files specified in a list
#'
#' @usage readQ.list(x = list("/path/to/QuantStudio1.xls", "/path/to/QuantStudio2.xls", "/path/to/QuantStudio3.xls"))
#' @param x A list of paths to a QuantStudio output file
#' @return A table of raw qPCR data from all imported QuantStudio output files
#'
#' @export
readQ.list <- function(filelist){

  outtable <- lapply(1:length(filelist), function(filenumber){

    data.frame(plate = filenumber, readQ.file(filelist[[filenumber]]))

  })

  tablenames <- table(unlist(lapply(outtable, names)))
  excludecols <- names(tablenames)[which(tablenames != length(outtable))]

  outtable <- lapply(outtable, function(y){

    ifelse(any(excludecols %in% names(y)),
           filteredtable <- y[,-(which(names(y) %in% excludecols))],
           filteredtable <- y)

    return(filteredtable)

    })

    outtable <- do.call(rbind, outtable)

    return(outtable)

}

#' Import QuantStudio output files
#'
#' This function imports raw qPCR data from a single QuantStudio output file path,
#' a list of QuantStudio output file paths, or a path to a folder containing QuantStudio output files
#'
#' @usage
#' readQ(x = "/path/to/QuantStudio.xls")
#' readQ(x = list("/path/to/QuantStudio1.xls", "/path/to/QuantStudio2.xls", "/path/to/QuantStudio3.xls"))
#' readQ(x = "/path/to/folder/")
#'
#' @param x Any of the following options: (1) a path to a single QuantStudio output file, (2) a list of paths to a QuantStudio output file, or (3) a path to a folder containing QuantStudio output files to be analyzed.
#' @return A table of raw qPCR data from all imported QuantStudio output files
#'
#' @export
readQ <- function(x){

  if(file_test("-d", x)){

    filelist <- list.files(x, pattern = "\\.xls*$")

    filelist <- lapply(filelist, function(filename){file.path(x, filename)})

    outtable <- readQ.list(filelist)

  } else if(is.list(x)){

    outtable <- readQ.list(x)

  } else if(file_test("-f", x)){

    outtable <- data.frame(plate = 1, readQ.file(x))

  }

  return(outtable)

}


#' Crosscheck Sample Names and Gene Names in the Uploaded Data Table
#'
#' This function produces a table of Sample Names and Gene Names from the uploaded
#' data table with a count of samples for each.
#'
#' @usage importCheck(datatable)
#'
#' @param datatable Imported data table from `readQ()` and associated functions
#' @return A table of sample names, gene names, and number of measurements for each
#'
#' @export
importCheck <- function(datatable){
  samplecount <- table(datatable[,c("Sample.Name", "Target.Name")])
  return(samplecount)
}
