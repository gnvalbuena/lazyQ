readQ.file <- function(x){

  filebase <- basename(x)

  outtable <- data.frame(file = filebase,
                         read_excel(path = x, sheet = "Results", skip = 46))

  outtable$CT[outtable$CT == "Undetermined"] <- NA
  outtable$CT <- as.numeric(outtable$CT)

  return(outtable)

}

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

