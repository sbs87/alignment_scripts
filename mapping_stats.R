#### ////////////////////////////////////////////////////////////
#### ////////////////////////////////////////////////////////////
#### mapping_stats.R
#### Author: Steven Smith
#### Version: 1
#### Date: 8/17/15
#### ------------------------------------------------------------
#### This script will compute basic post alignment mapping statistics for the SRL pipeline- includes mapping stats for rRNA, tRNA, etc & rarefaction. Also generates single counts table. 
#### ------------------------------------------------------------
#### Arguments:
#### 1- <file name for sample manifest>
#### 2- <the "n" used as the max # of mismatches to reference during BWA alignment
#### 3- <mirna count threshold- min read counts to call "feature detected">
#### 4- <subset feature- the feature NAME to match for the "subset feature" in rarefaction calculations
#### 5- <root_dir_path>
#### ------------------------------------------------------------
#### File/Dir Input requirements:
#### - trimmomatic trimming stats (formatted correctly)
#### - Sample name manifest
#### -  annotation files
#### -  alignment count files
#### ------------------------------------------------------------
#### Assumptions:
#### - Various input files have a particular file name format and location relative to root. This dependds on what the SRL pipeline outputs. 
#### - The reference sets are included in file names and also hardcoded.
#### - The piepine has summarized MIR and other SRL alignment counts before being input- this script simply adds it to the counts table. No fancy annotation footwork here. 
#### - The pipeline has summarized alignment stats and like the counts just adds it to a table.
#### ////////////////////////////////////////////////////////////
#### ////////////////////////////////////////////////////////////

## This version of the script exists on the grid- should be used for that envionment. 

## ////////////////////////////////////////////////////////////
## Setup enviornemnt
## ////////////////////////////////////////////////////////////

rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
#args <-c("SRL_manifest.txt",1,25,"MIR","/local/scratch2/steve/SRL/")
sample_names_fn<-args[1]#i.e., "SRL_sample_manifest.txt"
mismatch<-args[2]# 1
number_mirna_threshold<-as.numeric(args[3]) ## i.e, 25 counts- MAKE SURE NUMERIC!
subset_feature<-args[4] ##"MIR"
root_dir<-args[5] ## "/Users/stsmith/tmp/"
setwd(root_dir)

## ////////////////////////////////////////////////////////////
## Read in trimming stats and sample name manifest. Continue setting up tables
## ////////////////////////////////////////////////////////////

trimming_stats_fn<-paste0(root_dir,"SRL_trimmomatic.stats")
trimming_stats<-read.table(trimming_stats_fn,header=F,sep="\t")
names(trimming_stats)<-c("sample", "reads_in","reads_surviving","IGNORE_DROPPED","IRGNORE_DROPPED_PERCENT")
trimming_stats<-trimming_stats[,c("sample", "reads_in","reads_surviving")]

sample_names<-read.table(sample_names_fn,header=F,sep="\t")
sample_names<-sort(sample_names$V1)

## ////////////////////////////////////////////////////////////
## Setup mapping stats table
## ////////////////////////////////////////////////////////////

alignments<-c("tRNA","hum5SrDNA","HumRibosomal","Gvag","hg19") ## hardcoded- should be laoded
mapping_stats<-data.frame(matrix(0,nrow=length(sample_names),ncol=length(alignments)*2))
row.names(mapping_stats)<-sample_names
names(mapping_stats)<-paste(rep(alignments,each=2),rep(c("mapped","unmapped"),times=length(alignments)),sep=".")

## The following is a function that computes # features detected as a function of relative coverage. 
featureDetection_coverage<-function(relCoverage,sampleCount){
  sum(relCoverage*sampleCount>number_mirna_threshold)
}

## ////////////////////////////////////////////////////////////
## Setup rarefaction table
## ////////////////////////////////////////////////////////////

coverage<-seq(from = 1,to=0,by = -0.05)
rarefaction<-data.frame(matrix(0,nrow=length(coverage),ncol = 2*length(sample_names)))
row.names(rarefaction)<-coverage
names(rarefaction)<-paste(rep(c("global","subset"),times=length(sample_names)),rep(sample_names,each=2),sep=".")

## ////////////////////////////////////////////////////////////
## Setup counts table
## ////////////////////////////////////////////////////////////

counts_table<-data.frame(Feature="")


## ////////////////////////////////////////////////////////////
## Read in annotation counts and alignment stats tables to compute rarefaction, assign counts, and mapping stats. 
## ////////////////////////////////////////////////////////////

## For each sample, compute rarefaction at various relative coverages given the # of counts passing threshold 
for(i in 1:length(sample_names)){

  sample_name<-as.character(sample_names[i])
  ## Prepare counts table for sample for eventual merge to master counts table
  counts_path<-paste0(root_dir,"/annotation/SRL_",sample_name,"_R1_trimmed_gt16.hg19_n",mismatch,".hg19_annotation_miRNA_counts.txt")
  counts_table_sample<-read.table(counts_path,header=F,sep="\t")
  names(counts_table_sample)<-c(sample_name,"Feature")
  
  ## Compute rarefaction for all features in annotated SRL and also MIR (or "subset feature") only counts. 
  counts_table<-merge(counts_table,counts_table_sample,by.x = "Feature",by.y="Feature",all=T)
  counts_table_sample_subset<-counts_table_sample[grepl(subset_feature,counts_table_sample$Feature),]
  rarefaction[,paste0("global.",sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=counts_table_sample[1])    
  rarefaction[,paste0("subset.",sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=counts_table_sample_subset[1])
  
  ## Also store mapping statistics for each reference, i.e., for tRNA, rRNA, Gvag, etc. 
  for(j in 1:length(alignments)){
    alignment<-alignments[j]
    sample_path<-paste0(root_dir,"/alignment/",alignment,"_n",mismatch,"/SRL_",sample_name,"_R1_trimmed_gt16.aln_",alignment,"_n",mismatch,".stats.txt")
    sample_stats<-read.table(sample_path,header=F,sep="\t")
    aligned<-sum(sample_stats[!sample_stats$V1=="*","V3"])
    unaligned<-sample_stats[sample_stats$V1=="*","V4"]
    slot<-c(paste(alignment,"mapped",sep="."),paste(alignment,"unmapped",sep="."))
    mapping_stats[row.names(mapping_stats)==sample_name,slot]<-c(aligned,unaligned)
  }
}

## Move row names to a seperate column for easier merging to trimming stats
mapping_stats$sample<-row.names(mapping_stats)
mapping_stats<-merge(trimming_stats,mapping_stats,all=T) ## combine the reads that were trimmed from timmomatic to the mapping stats

## ////////////////////////////////////////////////////////////
## Summary stats for tables 
## ////////////////////////////////////////////////////////////

## Compute sumamry stats & percentages across all samples for mapping status. 
attach(mapping_stats)
total_input_tRNA<-tRNA.mapped + tRNA.unmapped
mapping_stats$total_input_tRNA<-total_input_tRNA
total_input_trimming<-reads_surviving
mapping_stats$unmapped_percent<-100*hg19.unmapped / (total_input_trimming)
mapping_stats$hg.mapped_ofInput_percent<-100*hg19.mapped / (total_input_trimming)
mapping_stats$input_difference<-total_input_trimming-total_input_tRNA
mapping_stats$mapped_total<-tRNA.mapped+hum5SrDNA.mapped+HumRibosomal.mapped+Gvag.mapped+hg19.mapped
detach(mapping_stats)

## Compute summary stats for rarefaction table
# for(s in c("global","subset")){
#   rarefaction_tmp<-rarefaction[,grepl(s,names(rarefaction))]
# rarefaction_mean<-apply(rarefaction,1,mean)
# rarefaction_stdev<-(apply(rarefaction,1,var))^(1/2)
# rarefaction_max<-apply(rarefaction,1,max)
# rarefaction_min<-apply(rarefaction,1,min)
# rarefaction_metrics<-data.frame(coverage=coverage,mean=rarefaction_mean,stdev=rarefaction_stdev,max=rarefaction_max,min=rarefaction_min)
# rarefaction_tmp$coverage<-coverage
# }

## ////////////////////////////////////////////////////////////
## Output tables to disk
## ////////////////////////////////////////////////////////////


rarefaction[is.na(rarefaction)]<-0 ##Cleanup rarefaction table
write.table(rarefaction,file=paste0(root_dir,"/stats/SRL_rarefaction_counts.gt",number_mirna_threshold,"_global_and_",subset_feature,"_n",mismatch,".txt"),sep="\t",row.names=T,col.names=NA,quote=F)
write.table(mapping_stats,file=paste0(root_dir,"/stats/SRL_mapping_statistics_n",mismatch,".txt"),sep="\t",row.names=F,quote=F)
write.table(counts_table,file=paste0(root_dir,"/stats/SRL_counts_table_n",mismatch,".txt"),sep="\t",row.names=F,quote=F)

## ////////////////////////////////////////////////////////////
## END OF SCRIPT
## ////////////////////////////////////////////////////////////

