suppressMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(snowfall)
  library(patchwork)
  source('code/utils.R')
  source('code/estimators.R')
})

files <- list.files(path = "code/experiments/configs", pattern = "*.R$", full.names = TRUE)

for( file in files ){
  print(file)
  plot_name <- paste("code/experiments/results/", gsub("\\.R$", "", basename(file)), ".png", sep="")
  if(!file.exists(plot_name)){
    if( grepl("mnar", basename(file)) ){
      source(file)
      source("code/experiments/simulation_mnar.R")
    }else if(grepl("questions", basename(file))){
      source(file)
      source("code/experiments/simulation_questions.R")
    }
  }
  
}
