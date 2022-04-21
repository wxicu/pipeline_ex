#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("utils"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
#suppressPackageStartupMessages(library("ggpubr"))

readDemux <- function(demux){
  demux <- fread(demux)
  demux_df <- demux[,c(2,5,6,7)]
  demux_df$cell.type  <- sapply(demux_df$BEST.GUESS,function(x){
    splitlist = strsplit(x,",")[[1]]
    if (splitlist[[1]] == splitlist[[2]]){ 
      splitlist[[2]]}
    else{
      "NA"}
  })
  demux_df[cell.type=="NA",]$cell.type = demux_df[cell.type=="NA",]$DROPLET.TYPE
  return (demux_df)
}

num_stat <- function(demux_df){
 #number statistics
  num_drop = nrow(demux_df)
  print(paste("number of droplets: ", num_drop ))
  num_sin = nrow(demux_df[demux_df$DROPLET.TYPE=="SNG",])
  prop_sin = round(num_sin*100/num_drop,2)
  print(paste("number of singlets: ", num_sin, ", with proportion ", prop_sin ,"%"))
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
# specify desired options, here take multiple result as input

#parser$add_argument("-t", "--tool", default="demuxlet", type="String", help="demultiplexing tool [default %(default)s]")
parser$add_argument("file", nargs='+', help="demultiplexing output file")

args <- parser$parse_args()
file <- args$file
#alpha <- args$alpha
#genoerror <-args$genoerror
result_df = as.data.frame(matrix(nrow= length(file), ncol=4))
colnames(result_df) = c("paramter","singlets","doublets","ambigous")

f1 = readDemux(file[1])
celltype = unique(append(f1$cell.type,c("AMB","DBL")))
num_drop = nrow(f1)
stat_celltype = as.data.frame(matrix(nrow= length(file), ncol=length(celltype)+1))
colnames(stat_celltype) = append("paramter",celltype)

assignment = as.data.frame(matrix(ncol= 1, nrow=num_drop))
assignment[,1] = f1$BARCODE
colnames(assignment) = c("barcode")
rindex = 1
for (f in file){
    demux_df = readDemux(f)
    num_drop = nrow(demux_df)
    num_sin = nrow(demux_df[demux_df$DROPLET.TYPE=="SNG",])
    prop_sin = round(num_sin*100/num_drop,2)
    num_dou = nrow(demux_df[demux_df$DROPLET.TYPE=="DBL",])
    prop_dou = round(num_dou*100/num_drop,2)
    num_amb = nrow(demux_df[demux_df$DROPLET.TYPE=="AMB",])
    prop_amb = round(num_amb*100/num_drop,2)
    row_name_element = strsplit(substr(f,10,nchar(f)-5), '_')[[1]]
    row_name = paste("alpha=",row_name_element[1],",error=",row_name_element[2],sep = "")
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb

    celltype = unique(demux_df$cell.type)
    for(c in celltype){
      num_cell = nrow(demux_df[demux_df$cell.type==c,])
      prop_cell =  round(num_cell/num_drop,2)
      stat_celltype[rindex,c] = prop_cell
    }
    stat_celltype[rindex,1] = row_name

    assignment[,row_name] = demux_df[order(assignment$barcode),]$DROPLET.TYPE
    rindex = rindex + 1
}
stat_celltype[is.na(stat_celltype)] = 0


write.table(result_df,"result.tsv",row.names=FALSE,sep="\t", quote = FALSE)
write.table(stat_celltype,"result_detail.tsv",row.names=FALSE,sep="\t", quote = FALSE)
write.table(assignment,"assignment.tsv",row.names=FALSE,sep="\t", quote = FALSE)

melt_df <- function(result){
  result_melt = melt(result)
  result_melt$paramter = sapply(result_melt$paramter,function(x){
    x = sub(',','\n',x)
    x = sub('alpha','a',x)
    x = sub('error','e',x)
  })
  return (result_melt)
}

bar_plot <- function(result_melt){
  barplot = ggplot(result_melt, aes( x = paramter, weight = value, fill = variable)) + 
    geom_bar(width = 0.7, position = "stack") + scale_fill_brewer(palette="Blues") + 
    theme_minimal() +ylab("Proportion")
  png("barplot.png")
  print(barplot)
  while (!is.null(dev.list()))  dev.off()
}

bar_plot_group <- function(result_melt){
  bar_plot_group = ggplot(result_melt, aes( x = paramter, y = value, fill = variable)) + 
    geom_bar(position="dodge",stat="identity") + 
    geom_text(aes(label=value), position=position_dodge(0.9),size = 4, vjust=-0.6) +
    scale_fill_brewer(palette="Blues")
    theme_minimal()
  png("barplot by group.png")
  print(bar_plot_group)
  while (!is.null(dev.list()))  dev.off() 
}

result_melt = melt_df(result_df)
result_detail_melt = melt_df(stat_celltype)
bar_plot(result_detail_melt)
bar_plot_group(result_melt)

