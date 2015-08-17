## This version of the script exists on the grid- should be used for that envionment. 

rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
#setwd("/local/projects-t2/HRBV/SRL")
setwd("/local/scratch2/steve/SRL")
root_dir<-getwd()
sample_names_fn<-args[1]#"SRL_manifest.sub.txt"
sample_names<-read.table(sample_names_fn,header=F,sep="\t")
sample_names<-sort(sample_names$V1)
#save(sample_names,file="SRL_sample_names.R")
#load("SRL_sample_names.R")
mismatch<-args[2]#0

alignments<-c("tRNA","hum5SrDNA","HumRibosomal","Gvag","hg19")
mapping_stats<-data.frame(matrix(0,nrow=length(sample_names),ncol=length(alignments)*2))
row.names(mapping_stats)<-sample_names
names(mapping_stats)<-paste(rep(alignments,each=2),rep(c("mapped","unmapped"),times=length(alignments)),sep=".")


#universal_set<-read.table("universal_set.txt",header=F,sep="\t")
#names(universal_set)<-"Feature"

#mir_set<-universal_set[grepl("MIR",universal_set$Feature),]

counts_table<-data.frame(Feature="")
number_mirna_threshold<-25

for(i in 1:length(sample_names)){
  #i<-2 
  sample_name<-as.character(sample_names[i])
  
  counts_path<-paste0("SRL_",sample_name,"_R1_trimmed_gt16.hg19_n",mismatch,".",alignment,"hg19_annotation_miRNA_counts.txt")
  counts_table_sample<-read.table(counts_path,header=F,sep="\t")
  names(counts_table_sample)<-c(sample_name,"Feature")
  
  counts_table<-merge(counts_table,new_table,by.x = "Feature",by.y="Feature")
  
  #unique_mirna_annotation_path<-paste0("~/Dropbox (IGS)/Working/smallRNA/Analysis/SRL_",sample_name,"_R1_trimmed_gt16.fastq.hg19.annotation.txt.uniqueMIRNAS.txt")
  #sample_rarefaction<-read.table(unique_mirna_annotation_path,header=F,sep="\t")
  ##rarefaction[,c(sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=new_table[1])
  #rarefaction_miRNA[,c(sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=new_table_mirna[1])
  
  
  for(j in 1:length(alignments)){
    #j<-3
    alignment<-alignments[j]
    sample_path<-paste0("alignment/",alignment,"_n",mismatch,"/SRL_",sample_name,"_R1_trimmed_gt16.aln_",alignment,"_n",mismatch,".stats.txt")
    sample_stats<-read.table(sample_path,header=F,sep="\t")
    aligned<-sum(sample_stats[!sample_stats$V1=="*","V3"])
    unaligned<-sample_stats[sample_stats$V1=="*","V4"]
    slot<-c(paste(alignment,"mapped",sep="."),paste(alignment,"unmapped",sep="."))
    mapping_stats[row.names(mapping_stats)==sample_name,slot]<-c(aligned,unaligned)
  }
}
mapping_stats$sample<-row.names(mapping_stats)


trimming_stats<-read.table("SRL_trimmomatic.stats",header=F,sep="\t")
names(trimming_stats)<-c("sample", "reads_in","reads_surviving","IGNORE_DROPPED","IRGNORE_DROPPED_PERCENT")
trimming_stats<-trimming_stats[,c("sample", "reads_in","reads_surviving")]

mapping_stats<-merge(trimming_stats,mapping_stats,all=T)

attach(mapping_stats)
total_input_tRNA<-tRNA.mapped + tRNA.unmapped
mapping_stats$total_input_tRNA<-total_input_tRNA
total_input_trimming<-reads_surviving
mapping_stats$unmapped_percent<-100*hg19.unmapped / (total_input_trimming)
mapping_stats$hg.mapped_ofInput_percent<-100*hg19.mapped / (total_input_trimming)
mapping_stats$input_difference<-total_input_trimming-total_input_tRNA
mapping_stats$mapped_total<-tRNA.mapped+hum5SrDNA.mapped+HumRibosomal.mapped+Gvag.mapped+hg19.mapped
detach(mapping_stats)

write.table(mapping_stats,file=paste0(root_dir,"/SRL_mapping_statistics_TEST_n",mismatch,".txt"),sep="\t",row.names=F,quote=F)
