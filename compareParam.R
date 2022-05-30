#!/usr/bin/Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("utils"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
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
  demux_df[cell.type=="NA",]$cell.type = "doublet"
  demux_df[DROPLET.TYPE=="AMB",]$cell.type = "ambigous"
  return (demux_df)
}


num_stat_demux <- function(demux_list){
  result_df = as.data.frame(matrix(nrow= length(demux_list), ncol=4))
  colnames(result_df) = c("Trial","singlets","doublets","ambigous/unassigned")
  
  f1 = readDemux(demux_list[1])
  
  celltype = unique(append(f1$cell.type,c("AMB","DBL")))
  stat_celltype = as.data.frame(matrix(nrow= length(demux_list), ncol=length(celltype)+1))
  colnames(stat_celltype) = append("Trial",celltype)
  
  num_drop = nrow(f1)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_drop))
  assignment[,1] = f1$BARCODE
  colnames(assignment) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=24, ncol=1))
  colnames(trial) = c("Parameter")
  trial[,1] = c("tag-group","tag-UMI","plp","field","geno-error-offset","geno-error-coeff","r2-info","min-mac","min-callrate","sm","sm-list","alpha","doublet-prior","sam-verbose","vcf-verbose","cap-BQ","min-BQ","min-MQ", "min-TD","excl-flag","group-list","min-total", "min-umi","min-snp")
  
  rindex = 1

  for (f in demux_list){
    demux_df = readDemux(f)
    num_sin = nrow(demux_df[demux_df$DROPLET.TYPE=="SNG",])
    prop_sin = round(num_sin*100/num_drop,2)
    num_dou = nrow(demux_df[demux_df$DROPLET.TYPE=="DBL",])
    prop_dou = round(num_dou*100/num_drop,2)
    num_amb = nrow(demux_df[demux_df$DROPLET.TYPE=="AMB",])
    prop_amb = round(num_amb*100/num_drop,2)
    
    row_name = paste("Demuxlet Trial",rindex)
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
    stat_celltype[is.na(stat_celltype)] = 0
    #question of order
    #assignment[,row_name] = demux_df[order(assignment$barcode),]$cell.type
    assignment[,row_name] = demux_df$cell.type
    
    #f = strsplit(f, split='/demux/',fixed=TRUE)
    f = strsplit(f, split='/demux+',fixed=TRUE)
    if(length(f[[1]]) == 1){
        f = f[[1]][1]
    }
    else{
	f = f[[1]][2]
    }
    #trial_param = strsplit(substr(f,7,nchar(f)-5), '+',fixed = TRUE)[[1]]
    trial_param = strsplit(f, '+',fixed = TRUE)[[1]]
    trial[,paste("Demuxlet Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df,"demuxlet_result.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(stat_celltype,"demuxlet_result_detail.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(assignment,"demuxlet_assignment.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(trial,"demuxlet_trial.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  
}


num_stat_vireo <-function(vireo_list){
  result_df = as.data.frame(matrix(nrow= length(vireo_list), ncol=4))
  colnames(result_df) = c("Trial","singlets","doublets","ambigous/unassigned")
  
  f1_summary = fread(paste(vireo_list[1],"/summary.tsv",sep = ""))
  f1_donorsid = fread(paste(vireo_list[1],"/donor_ids.tsv",sep = ""))
  
  donors = unique(append(f1_summary$Var1,c("doublet","unassigned")))
  stat_donors = as.data.frame(matrix(nrow= length(vireo_list), ncol=length(donors)+1))
  colnames(stat_donors) = append("Trial",donors)
  
  num_drop = sum(f1_summary$Freq)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_drop))
  
  assignment[,1] = f1_donorsid$cell
  colnames(assignment) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=15, ncol=1))
  colnames(trial) = c("Paramter")
  trial[,1] = c("n-donor","vartrix-data","donor-file","geno-tag","no-doublet","n-init","extra-donor",
                       "extra-donor-mode","force-learn-GT","ase-mode","no-plot","random-seed",
                       "cell-range","call-ambient-rna","n-proc")
  
  rindex = 1
  
  for (f in vireo_list){
    summary = fread(paste(f,"/summary.tsv",sep = ""))
    donorsid = fread(paste(f,"/donor_ids.tsv",sep = ""))
    
    num_dou = summary[summary$Var1 == 'doublet',]$Freq[1]
    prop_dou = round(num_dou*100/num_drop,4)
    num_amb = summary[summary$Var1 == 'unassigned',]$Freq[1]
    prop_amb = round(num_amb*100/num_drop,4)
    num_sin = num_drop - num_dou - num_amb
    prop_sin = round(num_sin*100/num_drop,4)
    
    row_name = paste("Vireo",rindex)
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    for(d in unique(summary$Var1)){
      num_cell = summary[summary$Var1 == d,]$Freq[1]
      prop_cell =  round(num_cell/num_drop,4)
      stat_donors[rindex,d] = num_cell
    }
    stat_donors[rindex,1] = row_name
    stat_donors[is.na(stat_donors)] = 0
    #assignment[,row_name] = donorsid[order(assignment$barcode),]$donor_id
    assignment[,row_name] = donorsid$donor_id
    
    f = strsplit(f, split='/vireo+',fixed=TRUE)
    #f = strsplit(f, split='/vireo/',fixed=TRUE)
    if(length(f[[1]]) == 1){
        f = f[[1]][1]
    }
    else{
	f = f[[1]][2]
    }
    #trial_param = strsplit(substr(f,7,nchar(f)), '+',fixed = TRUE)[[1]]
    trial_param = strsplit(f, '+',fixed = TRUE)[[1]]
    trial[, paste("Vireo Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
    
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df,"vireo_result.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(stat_donors,"vireo_result_detail.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(assignment,"vireo_assignment.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(trial,"vireo_trial.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  
  
}

num_stat_souporcell <-function(souporcell_list){
  result_df = as.data.frame(matrix(nrow= length(souporcell_list), ncol=4))
  colnames(result_df) = c("Trial","singlets","doublets","ambigous/unassigned")
  
  f1_cluster = fread(paste(souporcell_list[1],"/clusters.tsv",sep = ""))
  f1_cluster[f1_cluster$status == "doublet",]$assignment = "doublet"
  
  #f1_ambientRNA = fread(paste(souporcell_list[1],"/ambient_rna.txt",sep = ""))
  
  donors = unique(append(f1_cluster$assignment,c("doublet","unassigned")))
  
  stat_donors = as.data.frame(matrix(nrow= length(souporcell_list), ncol=length(donors)+1))
  colnames(stat_donors) = append("Trial",donors)
  
  num_drop = nrow(f1_cluster)
  assignments = as.data.frame(matrix(ncol= 1, nrow=num_drop))
  
  assignments[,1] = f1_cluster$barcode
  colnames(assignments) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=12, ncol=1))
  colnames(trial) = c("Parameter")
  trial[,1] = c("threads","clusters","ploidy","min-alt","min-ref","max-loci","restarts","with-common-variant","with-known-genotype","with-knowngenotypes-sample","skip-remap","ignore")
  
  rindex = 1
  
  for (f in souporcell_list){
    cluster = fread(paste(f,"/clusters.tsv",sep = ""))
    cluster[cluster$status == "doublet",]$assignment = "doublet"
    #ambientRNA = fread(paste(f,"/ambient_rna.txt",sep = ""))
    
    num_sin = nrow(cluster[cluster$status=="singlet",])
    prop_sin = round(num_sin*100/num_drop,2)
    num_dou = nrow(cluster[cluster$status=="doublet",])
    prop_dou = round(num_dou*100/num_drop,2)
    num_amb = nrow(cluster[cluster$status=="unassigned",])
    prop_amb = round(num_amb*100/num_drop,2)
    
    row_name = paste("Souporcell Trial",rindex)
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    donorlist = unique(append(cluster$assignment,c("doublet","unassigned")))
    
    for(d in donorlist){
      num_cell = nrow(cluster[cluster$assignment == d,])
      prop_cell =  round(num_cell/num_drop,4)
      stat_donors[rindex,d] = num_cell
    }
    stat_donors[rindex,1] = row_name
    stat_donors[is.na(stat_donors)] = 0
    
    #assignments[,row_name] = cluster[order(assignments$barcode),]$assignment
    assignments[,row_name] = cluster$assignment
    
    #f = strsplit(f, split='/souporcell/',fixed=TRUE)
    f = strsplit(f, split='/soup+',fixed=TRUE)
    if(length(f[[1]]) == 1){
        f = f[[1]][1]
    }
    else{
	f = f[[1]][2]
    }
    
    #trial_param = strsplit(substr(f,6,nchar(f)), '+',fixed = TRUE)[[1]]
    trial_param = strsplit(f, '+',fixed = TRUE)[[1]]
    trial[, paste("Souporcell Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
    
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df,"souporcell_result.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(stat_donors,"souporcell_result_detail.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(assignments,"souporcell_assignment.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(trial,"souporcell_trial.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  
  
}

num_stat_scSplit <-function(scSplit_list){
  result_df = as.data.frame(matrix(nrow= length(scSplit_list), ncol=4))
  colnames(result_df) = c("Trial","singlets","doublets","ambigous/unassigned")
  
  f1_result = fread(paste(scSplit_list[1],"/scSplit_result.csv",sep = ""))
  f1_prob = fread(paste(scSplit_list[1],"/scSplit_P_s_c.csv",sep = ""))
  
  donors = unique(append(f1_result$Cluster,c("unassigned")))
  stat_donors = as.data.frame(matrix(nrow= length(scSplit_list), ncol=length(donors)+1))
  colnames(stat_donors) = append("Trial",donors)
  
  num_drop = nrow(f1_result)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_drop))
  
  assignment[,1] = f1_result$Barcode
  colnames(assignment) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=6, ncol=1))
  colnames(trial) = c("Paramter")
  trial[,1] = c("with-common-snp","expected-num-samples","max-num-subpopulations-autodetect","num-EM","correction-doublets","known-genotypes")
  
  
  rindex = 1
  
  for (f in scSplit_list){
    result = fread(paste(f,"/scSplit_result.csv",sep = ""))
    prob = fread(paste(f,"/scSplit_P_s_c.csv",sep = ""))
    
    num_dou = nrow(result[substr(result$Cluster, 1, 3) %in% c("DBL"),])
    prop_dou = round(num_dou*100/num_drop,4)
    num_amb = 0
    prop_amb = 0
    num_sin = nrow(result[substr(result$Cluster, 1, 3) %in% c("SNG"),])
    prop_sin = round(num_sin*100/num_drop,4)
    
    row_name = paste("scSplit",rindex)
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    for(d in unique(result$Cluster)){
      num_cell = nrow(result[result$Cluster == d,])
      prop_cell =  round(num_cell/num_drop,4)
      stat_donors[rindex,d] = num_cell
    }
    stat_donors[rindex,1] = row_name
    stat_donors[is.na(stat_donors)] = 0
    
    #assignment[,row_name] = donorsid[order(assignment$barcode),]$donor_id
    assignment[,row_name] = result$Cluster
    
    f = strsplit(f, split='/scSplit+',fixed=TRUE)
    if(length(f[[1]]) == 1){
      f = f[[1]][1]
    }else{
      f = f[[1]][2]
    }
    trial_param = strsplit(f, '+',fixed = TRUE)[[1]]
    trial[, paste("scSplit Trial",rindex)] = trial_param 
  
    rindex = rindex + 1
    
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df,"scSplit_result.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(stat_donors,"scSplit_result_detail.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(assignment,"scSplit_assignment.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  write.table(trial,"scSplit_trial.tsv",row.names=FALSE,sep="\t", quote = FALSE)
  
  
}
#create parser object
parser <- ArgumentParser()
# specify desired options, here take multiple result as input

parser$add_argument("--tool", default="demuxlet", help="demultiplexing tool [default %(default)s]")
parser$add_argument("--file", nargs=1, help="demultiplexing output file")

args <- parser$parse_args()
file <- args$file
tool <- args$tool

demux_list = list()
vireo_list = list()
souporcell_list = list()
scSplit_list = list()

if (tool != "none"){
  file <- strsplit(file, ':',fixed = TRUE)[[1]]
  if (tool == "demuxlet"){
    demux_list= str_subset(file,'demux')
    print(demux_list)
    num_stat_demux(demux_list)
    
  }
  else if (tool == "vireo"){
    vireo_list = str_subset(file,'vireo')
    num_stat_vireo(vireo_list)
    
  }
  else if (tool == "souporcell"){
    souporcell_list = str_subset(file,'soup')
    num_stat_souporcell(souporcell_list)
    
  }
  else if (tool == "scSplit"){
    scSplit_list = str_subset(file,'scSplit')
    num_stat_scSplit(scSplit_list)
    
}
}else{
  demux_dir = paste(file,"/demux",sep="")
  demux_list = list.files (path = demux_dir, pattern="demux+", full.names=TRUE)
  vireo_dir = paste(file,"/vireo",sep="")
  vireo_list = list.files (path =  vireo_dir, pattern="vireo+", full.names=TRUE)
  soup_dir = paste(file,"/souporcell",sep="")
  souporcell_list = list.files (path =  soup_dir, pattern="soup+", full.names=TRUE)
  scSplit_dir = paste(file,"/scSplit",sep="")
  scSplit_list = list.files (path = scSplit_dir, pattern="scSplit+", full.names=TRUE)
  if (length(demux_list) != 0){
    num_stat_demux(demux_list)
  }else{
    print("No demuxlet output found!")
  }
  
  if (length(vireo_list) != 0){
    num_stat_vireo(vireo_list)
  }else{
    print("No vireo output found!")
  }
  
  if (length(souporcell_list) != 0){
    num_stat_souporcell(souporcell_list)
  }else{
    print("No souporcell output found!")
  }
  
  if (length(scSplit_list) != 0){
    num_stat_scSplit(scSplit_list)
  }else{
    print("No scSplit output found!")
  }

  assignment_all <- list.files(path =  ".", pattern = "_assignment.tsv", full.names = TRUE) %>%
    lapply(fread) %>%
    bind_cols()  %>% 
    as.data.frame()
  colnames(assignment_all)[1] = "Barcode"
  assignment_all = assignment_all[,-which(startsWith(colnames(assignment_all),'barcode'))]
  write.table(assignment_all,"assignment_all.tsv",row.names=FALSE,sep="\t", quote = FALSE)
}

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

#result_melt = melt_df(result_df)
#result_detail_melt = melt_df(stat_celltype)
#bar_plot(result_detail_melt)
#bar_plot_group(result_melt)


#file = "/home/xic/work/50/def4287076277e51cb85e99283a99f/demux+False+False+no_plp+GT+0.1+0.0+R2+1+0.50+False+no_sm_list+0.2+0.50+1000000+10000+40+13+20+0+3844+False+0+0+0.best:/home/xic/work/89/d427f3e2b4b230fb18ab15c8c78ee9/demux+False+False+no_plp+GT+0.1+0.0+R2+1+0.50+False+no_sm_list+0.5+0.50+1000000+10000+40+13+20+0+3844+False+0+0+0.best"
#file = "/home/xic/work/83/ff58880f4a2354ac1030b5c84c9581/vireo+4+no_vartrixData+no_donorfile+PL+False+50+0+distance+False+False+False+none+all+False+1"

