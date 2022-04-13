library(R.utils)
library(ggplot2);
library(data.table);
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
#else if (length(args)==1) {
  # default output file
 # args[2] = "out.txt"
#}

demux <- fread(args[1]);

demux_df <- demux[,c(2,5,6,7)]
demux_df$cell.type  <- sapply(demux_df$BEST.GUESS,function(x){
  splitlist = strsplit(x,",")[[1]]
  if (splitlist[[1]] == splitlist[[2]]){ 
    splitlist[[2]]}
  else{
    "DBL"}
})

#violin plot of log likelihood
vlplot <- ggplot(demux_df, aes(cell.type, BEST.LLK, color = cell.type)) + geom_violin() + geom_boxplot(width=0.1, position = position_dodge(0.9)) + ylab("log likelihood")+ xlab("different cell types")
png("violinplot.png")
print(vlplot)
while (!is.null(dev.list()))  dev.off()


# number statistics
num_drop = nrow(demux_df)
print(paste("number of droplets: ", num_drop ))
num_sin = nrow(demux_df[demux_df$DROPLET.TYPE=="SNG",])
prop_sin = round(num_sin*100/num_drop,2)
print(paste("number of singlets: ", num_sin, ", with proportion ",prop_sin ,"%"))
num_dou = nrow(demux_df[demux_df$DROPLET.TYPE=="DBL",])
prop_dou = round(num_dou*100/num_drop,2)
print(paste("number of doublets: ", num_dou, ", with proportion ", prop_dou,"%"))
celltyple = unique(demux_df$cell.type)
for( c in celltyple ){
  if (c!="DBL"){
    num_cell = nrow(demux_df[demux_df$cell.type==c,])
    prop_cell =  round(num_cell*100/num_drop,2)
    print(paste("Among the singlets, number of ", c , ": ", num_cell, ", with proportion ", prop_cell,"%"))
    
  }
}

