library(R.utils)
library(ggplot2);
library(data.table);
args = commandArgs(trailingOnly=TRUE)

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
dev.off()


# number statistics
print(paste("number of singlets: ", nrow(demux_df[demux_df$DROPLET.TYPE=="SNG",])))
print(paste("number of doublets: ", nrow(demux_df[demux_df$DROPLET.TYPE=="DBL",])))
print("Among the singlets,")
celltyple = unique(demux_df$cell.type)
for( c in celltyple ){
  if (c != "DBL"){
	print(paste("number of ", c , ": ",nrow(demux_df[demux_df$cell.type==c,])))
}
  
}

