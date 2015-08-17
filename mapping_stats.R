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
subset_feature<-args[4]
alignments<-c("tRNA","hum5SrDNA","HumRibosomal","Gvag","hg19")
mapping_stats<-data.frame(matrix(0,nrow=length(sample_names),ncol=length(alignments)*2))
row.names(mapping_stats)<-sample_names
names(mapping_stats)<-paste(rep(alignments,each=2),rep(c("mapped","unmapped"),times=length(alignments)),sep=".")

featureDetection_coverage<-function(relCoverage,sampleCount){
  sum(relCoverage*sampleCount>number_mirna_threshold)
}

coverage<-seq(from = 1,to=0,by = -0.05)
rarefaction<-data.frame(matrix(0,nrow=length(coverage),ncol = 2*length(sample_names)))
row.names(rarefaction)<-coverage
names(rarefaction)<-paste(rep(c("global","subset"),times=length(sample_names)),rep(sample_names,each=2),sep=".")
rarefaction_miRNA<-rarefaction

#mir_set<-universal_set[grepl("MIR",universal_set$Feature),]

counts_table<-data.frame(Feature="")
number_mirna_threshold<-args[3]#25

for(i in 1:length(sample_names)){
  #i<-2
	
  sample_name<-as.character(sample_names[i])
 
    counts_path<-paste0(root_dir,"/annotation/SRL_",sample_name,"_R1_trimmed_gt16.hg19_n",mismatch,".hg19_annotation_miRNA_counts.txt")
    counts_table_sample<-read.table(counts_path,header=F,sep="\t")
    names(counts_table_sample)<-c(sample_name,"Feature")
    
    counts_table<-merge(counts_table,counts_table_sample,by.x = "Feature",by.y="Feature",all=T)
	counts_table_sample_subset<-counts_table_sample[grepl(subset_feature,counts_table_sample$Feature),]
rarefaction[,paste0("global.",sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=counts_table_sample[1])    
rarefaction[,paste0("subset.",sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=counts_table_sample_subset[1])
 
  for(j in 1:length(alignments)){
    #j<-3
    alignment<-alignments[j]
    sample_path<-paste0(root_dir,"/alignment/",alignment,"_n",mismatch,"/SRL_",sample_name,"_R1_trimmed_gt16.aln_",alignment,"_n",mismatch,".stats.txt")
    sample_stats<-read.table(sample_path,header=F,sep="\t")
    aligned<-sum(sample_stats[!sample_stats$V1=="*","V3"])
    unaligned<-sample_stats[sample_stats$V1=="*","V4"]
    slot<-c(paste(alignment,"mapped",sep="."),paste(alignment,"unmapped",sep="."))
    mapping_stats[row.names(mapping_stats)==sample_name,slot]<-c(aligned,unaligned)
    
    
  }
}
mapping_stats$sample<-row.names(mapping_stats)
rarefaction

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
