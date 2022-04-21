#!/usr/bin/env Rscript
#install.packages("argparse")
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))


readDemux <- function(demux){
  demux_df <- demux[,c(2,5,6,7)]
  demux_df$cell.type  <- sapply(demux_df$BEST.GUESS,function(x){
    splitlist = strsplit(x,",")[[1]]
    if (splitlist[[1]] == splitlist[[2]]){ 
      splitlist[[2]]}
    else{ "NA"}
  })
  
  demux_df[cell.type=="NA",]$cell.type = demux_df[cell.type=="NA",]$DROPLET.TYPE
  return (demux_df)
}


num_stat <- function(demux_df){
  # number statistics
  num_drop = nrow(demux_df)
  print(paste("number of droplets: ", num_drop ))
  num_sin = nrow(demux_df[demux_df$DROPLET.TYPE=="SNG",])
  prop_sin = round(num_sin*100/num_drop,2)
  print(paste("number of singlets: ", num_sin, ", with proportion ",prop_sin ,"%"))
  num_dou = nrow(demux_df[demux_df$DROPLET.TYPE=="DBL",])
  prop_dou = round(num_dou*100/num_drop,2)
  print(paste("number of doublets: ", num_dou, ", with proportion ", prop_dou,"%"))
  num_amb = nrow(demux_df[demux_df$DROPLET.TYPE=="AMB",])
  prop_amb = round(num_amb*100/num_drop,2)
  print(paste("number of ambigous droplets: ", num_amb, ", with proportion ", prop_amb,"%"))
  
  celltyple = unique(demux_df$cell.type)
  for( c in celltyple ){
    if (c!="DBL" && c!="AMB"){
      num_cell = nrow(demux_df[demux_df$cell.type==c,])
      prop_cell =  round(num_cell*100/num_drop,2)
      print(paste("Among the singlets, number of ", c , ": ", num_cell, ", with proportion ", prop_cell,"%"))
    }
  }
  
}



#create parser object
parser <- ArgumentParser()
# specify desired options 
#didnt decide whether to allow only one output file or multiple
#here just take one single result as input

parser$add_argument("-t", "--tool", default="demuxlet", help="demultiplexing tool [default %(default)s]")
parser$add_argument("file", nargs=1, help="demultiplexing output file")

args <- parser$parse_args()
file <- args$file 
if( file.access(file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file))
} 

if( args$tool == "demuxlet") {
  demux = fread(file)
  demux_df = readDemux(demux)
  num_stat(demux_df)
  
  #violin plot of log likelihood
  vlplot <- ggplot(demux_df, aes(cell.type, BEST.LLK, color = cell.type)) + geom_violin() + geom_boxplot(width=0.1, position = position_dodge(0.9)) + ylab("log likelihood")+ xlab("different cell types")
  png("violinplot_demux.png")
  print(vlplot)
  while (!is.null(dev.list()))  dev.off()
} 




