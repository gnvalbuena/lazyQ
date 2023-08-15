#' @export
analyzeQ.CT <- function(readqtable, housekeepinggene){

  workingtable <- readqtable

  samples <- unique(workingtable$Sample.Name)
  samples <- samples[order(samples)]
  targets <- unique(workingtable$Target.Name)
  targets <- targets[which(targets != housekeepinggene)]
  targets <- targets[order(targets)]

  n.housekeep <- max(table(workingtable[workingtable$Target.Name == housekeepinggene, c("Sample.Name", "Target.Name")]))
  n.targets <- max(table(workingtable[workingtable$Target.Name != housekeepinggene, c("Sample.Name", "Target.Name")]))

  # Prepare housekeeping gene table

  housekeep.table <- lapply(samples, function(samplename){

    sample.subtable <- workingtable[which(workingtable$Sample.Name == samplename & workingtable$Target.Name == housekeepinggene),
                                    c("Sample.Name", "Target.Name", "CT", "Ct.Mean", "Ct.SD")]

    outtable.colnames <- c("Sample.Name", "Housekeeping.Gene", paste0("housekCT.", 1:n.housekeep),
                           "housekCT.mean", "housekCT.sd")

    outtable <- data.frame(data.frame(matrix(nrow = 0, ncol = length(outtable.colnames))))
    colnames(outtable) <- outtable.colnames

    outtable[1,"Sample.Name"] <- samplename
    outtable[1, "Housekeeping.Gene"] <- housekeepinggene

    for(x in 1:nrow(sample.subtable)){
      outtable[1,(x+2)] <- sample.subtable$CT[x]
    }

    outtable$housekCT.mean <- mean(sample.subtable$CT, na.rm = TRUE)
    outtable$housekCT.sd <- sd(sample.subtable$CT, na.rm = TRUE)

    return(outtable)

  })

  housekeep.table <- do.call(rbind, housekeep.table)

  # Prepare target table and append to corresponding housekeeping gene values

  targettable <- lapply(targets, function(targetname){

    out.target <- lapply(as.list(housekeep.table$Sample.Name), function(samplename){

      sample.subtable <- workingtable[which(workingtable$Sample.Name == samplename & workingtable$Target.Name == targetname),
                                      c("Sample.Name", "Target.Name", "CT", "Ct.Mean", "Ct.SD")]

      outtable.colnames <- c("Sample.Name", "Target.Gene", paste0("targetCT.", 1:n.targets),
                             "targetCT.mean", "targetCT.sd")

      outtable <- data.frame(data.frame(matrix(nrow = 0, ncol = length(outtable.colnames))))
      colnames(outtable) <- outtable.colnames

      outtable[1,"Sample.Name"] <- samplename
      outtable[1, "Target.Gene"] <- targetname

      for(x in 1:nrow(sample.subtable)){
        outtable[1,(x+2)] <- sample.subtable$CT[x]
      }

      outtable$targetCT.mean <- mean(sample.subtable$CT, na.rm = TRUE)
      outtable$targetCT.sd <- sd(sample.subtable$CT, na.rm = TRUE)

      return(outtable)

    })

    out.target <- do.call(rbind, out.target)

    out.target <- cbind(housekeep.table, out.target)

  })

  targettable <- do.call(rbind, targettable)
  targettable$dCT <- targettable$targetCT.mean - targettable$housekCT.mean

  return(targettable)

}

#' @export
analyzeQ.ddCT <- function(target.table, control.samples){

  workingtable <- cbind(Sample.Name = target.table[,1],
                        sample.class = ifelse(target.table$Sample.Name %in% control.samples, "Control", NA),
                        target.table[,2:ncol(target.table)])

  workingtable <- do.call(rbind, lapply(unique(workingtable$Target.Gene), function(targetgene){

    target.subtable <- subset(workingtable, workingtable$Target.Gene == targetgene)
    target.subtable$control.dCT.mean <- mean(target.subtable$dCT[target.subtable$sample.class == "Control"], na.rm = TRUE)
    target.subtable$control.dCT.sd <- sd(target.subtable$dCT[target.subtable$sample.class == "Control"], na.rm = TRUE)
    return(target.subtable)

  }))

  workingtable$ddCT <- workingtable$dCT - workingtable$control.dCT.mean
  workingtable$log2FC <- -(workingtable$ddCT)
  workingtable$fold.change <- 2^(-(workingtable$ddCT))

  return(workingtable)

}

#' @export
analyzeQ <- function(readqtable, housekeepinggene, control.samples){

  target.table <- analyzeQ.CT(readqtable, housekeepinggene)

  ddCT.table <- analyzeQ.ddCT(target.table, control.samples)

  return(ddCT.table)

}






