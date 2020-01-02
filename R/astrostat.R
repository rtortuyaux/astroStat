#Author: romain_tortuyaux (rtortuyaux@live.fr)
#GitHub: rtortuyaux/astroStat
#College de France, Paris, France
#07 juin 2019#

#####          AstroStat - astroDot data analysis          #####

#' astrostat
#'
#' @param a
#' @param b
#' @param paired
#' @param graph
#' @param ratio
#'
#' @return
#' @export
#'
#' @examples

astrostat <- function(a = "WT", b = "Ex", paired = FALSE, ratio = FALSE, graph = FALSE)
{

  #install missing libraries
  requiredPackages = c("base", "tools", "plyr", "dplyr", "ggplot2", "rstudioapi", "svDialogs")
  for(i in requiredPackages){
    if(!require(i,character.only = TRUE)) install.packages(i)
    library(i,character.only = TRUE)
  }

  #check sysname
  linux = FALSE
  if (Sys.info()["sysname"] == "Linux") {
    linux = TRUE
  }

  #template analysis
  print("Welcome in astroStat : astroDot data analysis")
  if (!linux) {
    print("Choose working directory")
    print("Press enter to continue")
    pause <- readLines (n=1)
    Dir <- selectDirectory(caption = "Select astroDot data directory", label = "Select", path = NULL)
    setwd(Dir)
    form <- list(
      "Condition1:TXT" = "WT",
      "Condition2:TXT" = "Ex",
      "Paired:CHK" = FALSE,
      "Ratio_analysis:CHK" = FALSE,
      "Normality_plot:CHK" = FALSE)
    parameters <- dlg_form(form, "astroStat parameters")$res
    a <- parameters$Condition1
    b <- parameters$Condition2
    paired <- as.logical(parameters$Paired)
    ratio <- as.logical(parameters$Ratio_analysis)
    graph <- as.logical(parameters$Normality_plot)
  } else {
    form <- list(
      "Directory:DIR" = getwd(),
      "Condition1:TXT" = "WT",
      "Condition2:TXT" = "Ex",
      "Paired:CHK" = FALSE,
      "Ratio_analysis:CHK" = FALSE,
      "Normality_plot:CHK" = FALSE)
    parameters <- dlg_form(form, "astroStat parameters")$res
    Dir <- parameters$Directory
    a <- parameters$Condition1
    b <- parameters$Condition2
    paired <- as.logical(parameters$Paired)
    ratio <- as.logical(parameters$Ratio_analysis)
    graph <- as.logical(parameters$Normality_plot)
    setwd(Dir)
  }

  #look for astroDots files in working directory
  if(length(list.files(pattern = "AstroDot")) != 0) {
    cat("astroDot files found", "\n")
    cat("astroStat in process...", "\n")
    Sys.sleep(1)
  } else {
    stop("astroDot files not found : end of processus")
  }

  #open and bind files according to group
  g <- c(a, b)
  files_names <- vector ("list", length = 2)
  data_names <- vector("list", length = length((files_names)))

  for (i in seq_along(g)) {
    files_names[[i]] <- list.files(path = getwd(), pattern = g[[i]])
    data_names[[i]] <- file_path_sans_ext(files_names[[i]])
    for (n in seq_along(data_names[[i]])) {
      assign(data_names[[i]][n], read.csv(files_names[[i]][n], sep = "\t"))
    }
  }

  #check summary of each mouse
  for (i in seq_along(g)) {
    for (j in 1 : length(data_names[[i]])) {
      print(data_names[[i]][j])
      print(summary(get(data_names[[i]][j])))
    }
  }

  data <- vector("list", length = length(g))
  data[[1]] <- get(data_names[[1]][1])
  data[[2]] <- get(data_names[[2]][1])

  #verification
  if (length(files_names[[1]]) != length(files_names[[2]]))
    stop("different file number between conditions - end of processus")

  if (length(data_names[[1]]) > 1) {
    for (i in seq_along(g)) {
      for (n in 2 : length(data_names[[1]])) {
        data[[i]] <- rbind(data[[i]], get(data_names[[i]][n]))
      }
    }
  }

  #verification
  if (paired && (nrow(data[[1]]) != nrow(data[[2]])))
    stop("different astrocyte number between conditions - end of processus")
  if (!paired && ratio) {
    ratio = FALSE
    cat("ratio comparison on unpaired analysis is not possible", "\n")
  }

  #create .txt result file
  wd <- paste0("AstroStat_", substr(basename(files_names[[1]][1]), 14, 18), substr(basename(files_names[[1]][1]), 19, 22), substr(basename(files_names[[1]][1]), 28, 30), "___", substr(basename(files_names[[2]][1]), 14, 18), substr(basename(files_names[[2]][1]), 19, 22), substr(basename(files_names[[2]][1]), 28, 30))
  dir.create(wd)
  setwd(wd)
  dir.create("NormaTest")
  sink(paste0(wd, ".txt"), split = TRUE)

  cat("###astroStat : astroDots data analysis###", "\n")
  cat("<author> romain_tortuyaux (rtortuyaux@live.fr)", "\n")
  cat("College de France, Paris, France", "\n", "\n")
  cat("Performed on :", date(), "by", Sys.info()['user'], "\n")
  cat("Working directory :", Dir, "\n")

  #Parameters of analysis
  cat("Compared conditions :", a, b, "\n")
  cat("Paired analysis :", paired, "\n")
  cat("normality plot :", graph, "\n", "\n")

  #loaded files
  print("loaded files")
  cat("\n")
  print(paste("group", a))
  cat("processed files :", basename(files_names[[1]]), "\n", "\n")
  print(paste("group", b))
  cat("processed files :", basename(files_names[[2]]), "\n")

  cat("\n", "\n")
  print("DESCRIPTIVE STATISTICS")
  cat("\n")

  for (i in seq_along(g)) {
    print(paste("group", g[[i]]))
    for (n in 1 : length(files_names[[i]])) {
      cat("Astrocyte number in mouse", substr(basename(files_names[[i]][n]), 24, 26), ":", nrow(get(data_names[[i]][n])), "\n")
    }
  }

  cat("\n", "\n")
  cat("Astrocyte number in group", a, ":", nrow(data[[1]]), "\n")
  cat("Astrocyte number in group", b, ":", nrow(data[[2]]), "\n")

  if (paired) {
    tableRatio <- data.frame(
      ratio = data[[1]]$Density.dots.in.Astro / data[[2]]$Density.dots.in.Astro
    )
    write.table(tableRatio, file = paste0(getwd(), "/tableRatio_", wd, ".csv"), row.names = FALSE, sep = "\t")
  }

  #Remove unused column
  speRNA = FALSE
  if (sum(data[[i]]$Percentage.of.dot.in.GFAP.neg != 0)) {
    speRNA = TRUE
  }
  for (i in seq_along(data)) {
    data[[i]]$ImageName <- NULL
    data[[i]]$Roi.Name <- NULL
    data[[i]]$bg.Intensity <- NULL
    if (paired) {
      data[[i]]$Astrocyte.Volume <- NULL
    }
    if (!speRNA) {
      data[[i]]$Percentage.of.dot.in.GFAP.neg <- NULL
    } else {
      data[[i]]$Density.dots.in.Astro <- NULL
    }
  }

    #descriptive statistics
  for (i in seq_along(g)) {
    cat("\n", "\n")
    print(paste("descriptive statistics of", g[i]))
    for (n in names(data[[i]])) {
      cat("\n", "\n", colnames(data[[i]][n]), "\n", "\n")
      cat(summary(data[[i]][n]), "\n", "\n")
      cat("sd =", sd(data[[i]][[n]]), "\n")
    }
  }

  #process mean test
  #print results in .csv table
  cat("\n", "\n")
  print("mean comparison")
  cat("data : x =", a, "; y =", b, "\n")
  cat("\n")
  tableRes <- data.frame()
  for (i in names(data[[1]])) {
    colname <- colnames(data[[1]][i])
    print(colname)
    res = meanTest(data[[1]][[i]], data[[2]][[i]], paired = paired, graph = graph, colname = colnames(data[[1]][i]))
    print(res)
    if ( i == "Density.dots.in.Astro") {
      data_i <- data.frame(
        X = paste0( round( mean(data[[1]][[i]]), 4) * 100, " (", round( sd(data[[1]][[i]]), 4) * 100, ")"),
        Y = paste0( round( mean(data[[2]][[i]]), 4) * 100, " (", round( sd(data[[2]][[i]]), 4) * 100, ")"),
        pvalue = res$p.value,
        row.names = paste(colname, "(*100)")
      )
      tableRes <- rbind(tableRes, data_i)
      } else {
      data_i <- data.frame(
        X = paste0(round(mean(data[[1]][[i]]), 2), " (", round(sd(data[[1]][[i]]), 2), ")"),
        Y = paste0(round(mean(data[[2]][[i]]), 2), " (", round(sd(data[[2]][[i]]), 2), ")"),
        pvalue = res$p.value,
        row.names = colname
      )
      tableRes <- rbind(tableRes, data_i)
    }
  }
  names(tableRes)[1] <- paste0(a, " (n = ", nrow(data[[1]]), ")")
  names(tableRes)[2] <- paste0(b, " (n = ", nrow(data[[2]]), ")")
  write.csv(tableRes, file = paste0(getwd(), "/tableRes_", wd, ".csv"), row.names = TRUE)
  print(tableRes)

  #make camembert
  dir.create("astroPlot")
  if (!paired) {
  for (i in seq_along(g)) {
  df <- data.frame(
    localisation = c("soma", "large branche", "fine branche"),
    value = c(round(mean(data[[i]][[4]]), 2), round(mean(data[[i]][[6]]), 2), round(mean(data[[i]][[5]]), 2))
  )
  bp <- ggplot(df, aes (x = "", y=value, fill=localisation)) +
          geom_bar(width = 1, stat = "identity")
  camembert <- bp + coord_polar("y", start = 0) +
    theme_bw(base_size = 16) +
    labs(x = element_blank(), y = element_blank()) +
    ggtitle(paste("Mean distribution of mRNA in astrocyte \n", g[i])) +
    theme(plot.title = element_text(hjust = 0.5, size=20, face = "bold", color = "blue"),
          axis.ticks = element_blank(),
          axis.text.x=element_blank(),
          panel.grid=element_blank(),
          legend.title = element_text(face="bold")) +
    geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]),
                  label = paste(round(value, 0), "%")), size=7)
  ggsave(paste0("astroCam_",g[i],".pdf"), plot = camembert, path = paste0(getwd(), "/astroPlot"),
         units="mm",width=250, height = 250, dpi=300)
  }
}

  # sum column (fine and large processes)
  # define new variable (process) + meantest
  cat("\n", "\n")
  print("comparison of dot in processes (fine + large)")
  process_WT <- data.frame(
    process = data[[1]]$Percentage.of.dots.in.fine.processes + data[[1]]$Percentage.of.dots.in.large.processes
  )
  process_APP <- data.frame(
    process = data[[2]]$Percentage.of.dots.in.fine.processes + data[[2]]$Percentage.of.dots.in.large.processes
  )
  print(paste("mean (sd) of", a))
  print(paste0(round(mean(process_WT$process), 2), " (", round(sd(process_WT$process), 2), ")"))
  print(paste("mean (sd) of", b))
  print(paste0(round(mean(process_APP$process), 2), " (", round(sd(process_APP$process), 2), ")"))

  print(meanTest(process_WT$process, process_APP$process, paired = FALSE, graph = FALSE))

  #create a .csv with all data
  for (i in seq_along(g))
    if (i == 1) {
      condition <- a
      data[[i]] <- cbind(data[[i]], condition)
    } else {
      condition <- b
      data[[i]] <- cbind(data[[i]], condition)
    }
  data_synthese <- rbind(data[[1]], data[[2]])
  write.table(data_synthese, file = paste0(getwd(), "/data_synthese_", wd, ".csv"), row.names = FALSE, sep = "\t")

  #create plot (ggplot2)
  condition <- as.factor(data_synthese$condition)
  for (i in colnames(data_synthese)) {
    if (i == "condition") {
      break
    }
    text_title <- gsub("\\."," ",i)
    text_x <- gsub("Astrocyte"," ",text_title)
    if (i == "Density.dots.in.Astro") {
      meanDistrib <- ddply(data_synthese, "condition", summarise, grp.mean=mean(Density.dots.in.Astro)*100)
      histo <- ggplot(data_synthese, aes(data_synthese[[i]]*100, color=condition, fill=condition)) +
        theme_bw(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, size=20, face = "bold", color = "blue"),
              axis.title.x = element_text(face="bold"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.title = element_text(face="bold")) +
        geom_histogram(aes(y= ..density..), position="identity", alpha=0.5, bins = 30) +
        geom_density(alpha=0.3, show.legend = FALSE) +
        labs(x = "mRNA density (*100)", y = "") +
        geom_vline(data=meanDistrib, aes(xintercept=grp.mean, color=condition), linetype="dashed", size = 2, show.legend = FALSE)
        ggtitle(paste("Distribution histogram \n --", text_title, "--"))
        } else {
      histo <- ggplot(data_synthese, aes(data_synthese[[i]], color=condition, fill=condition)) +
        theme_bw(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, size=20, face = "bold", color = "blue"),
              axis.title.x = element_text(face="bold"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.title = element_text(face="bold")) +
        geom_histogram(aes(y= ..density..), position="identity", alpha=0.5, bins = 30) +
        geom_density(alpha=0.3, show.legend = FALSE) +
        labs(x = text_title, y = "") +
        ggtitle(paste("Distribution histogram \n --", text_title, "--"))
        }
      ggsave(paste0("astroPlot_",i,".pdf"), plot = histo, path = paste0(getwd(), "/astroPlot"),
           units="mm",width=250, height = 250, dpi=300)
  }

  if(paired) {
    cat("\n", paste0("mean ratio (SD) of mRNA density ", a, " / ", b, " = ", round(mean(tableRatio[[1]]), 3), " (", round(sd(tableRatio[[1]]), 3), ")"), "\n", "\n")
  }

  #Comparison of ratio in paired analysis
  if(paired && ratio) {
    sink()
    print("Comparison of ratio : select .csv file")
    print("Press enter to continue")
    pause <- readLines (n=1)
    tableRatioC_F <- file.choose()
    tableRatioC <- read.csv(tableRatioC_F, sep = "\t")
    sink(paste0(wd, ".txt"), split = TRUE, append = TRUE)
    print("Comparison of ratio :")
    cat(wd, "versus", basename(tableRatioC_F), "\n", "\n")
    resRatio <- meanTest(tableRatio$ratio, tableRatioC$ratio, paired = FALSE, graph = graph, colname = colnames(tableRatioC))
    print(resRatio)
  }

  if(speRNA) {
    cat("\n", "exclusion ratio not processed due to astroDots parameters", "\n", "\n")
  }

  #clear console
  cat("\014", "\n")
  cat("/process completed.")
  sink()
  file.show(paste0(wd, ".txt"), title = paste0(wd, ".txt"))
  remove(list = ls())
}
