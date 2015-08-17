## This version of the script exists on the grid- should be used for that envionment. 

rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
sample_names_fn<-args[1]#i.e., "SRL_manifest.txt"
mismatch<-args[2]
number_mirna_threshold<-args[3] ## i.e, 25 counts
subset_feature<-args[4]
root_dir<-args[5]
setwd(root_dir)

trimming_stats_fn<-paste0(root_dir,"SRL_trimmomatic.stats")
trimming_stats<-read.table(trimming_stats_fn,header=F,sep="\t")
names(trimming_stats)<-c("sample", "reads_in","reads_surviving","IGNORE_DROPPED","IRGNORE_DROPPED_PERCENT")
trimming_stats<-trimming_stats[,c("sample", "reads_in","reads_surviving")]

sample_names<-read.table(sample_names_fn,header=F,sep="\t")
sample_names<-sort(sample_names$V1)

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

counts_table<-data.frame(Feature="")

for(i in 1:length(sample_names)){

  sample_name<-as.character(sample_names[i])
  
  counts_path<-paste0(root_dir,"/annotation/SRL_",sample_name,"_R1_trimmed_gt16.hg19_n",mismatch,".hg19_annotation_miRNA_counts.txt")
  counts_table_sample<-read.table(counts_path,header=F,sep="\t")
  names(counts_table_sample)<-c(sample_name,"Feature")
  
  counts_table<-merge(counts_table,counts_table_sample,by.x = "Feature",by.y="Feature",all=T)
  counts_table_sample_subset<-counts_table_sample[grepl(subset_feature,counts_table_sample$Feature),]
  rarefaction[,paste0("global.",sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=counts_table_sample[1])    
  rarefaction[,paste0("subset.",sample_name)]<-sapply(coverage,featureDetection_coverage,sampleCount=counts_table_sample_subset[1])
  
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

mapping_stats$sample<-row.names(mapping_stats)
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

# for(s in c("global","subset")){
#   rarefaction_tmp<-rarefaction[,grepl(s,names(rarefaction))]
# rarefaction_mean<-apply(rarefaction,1,mean)
# rarefaction_stdev<-(apply(rarefaction,1,var))^(1/2)
# rarefaction_max<-apply(rarefaction,1,max)
# rarefaction_min<-apply(rarefaction,1,min)
# rarefaction_metrics<-data.frame(coverage=coverage,mean=rarefaction_mean,stdev=rarefaction_stdev,max=rarefaction_max,min=rarefaction_min)
#rarefaction$coverage<-coverage
# }
rarefaction[is.na(rarefaction)]<-0
write.table(rarefaction,file=paste0(root_dir,"/stats/SRL_rarefaction_counts.gt",number_mirna_threshold,"_global_and_",subset_feature,"_n",mismatch,".txt"),sep="\t",row.names=T,col.names=NA,quote=F)

write.table(mapping_stats,file=paste0(root_dir,"/stats/SRL_mapping_statistics_n",mismatch,".txt"),sep="\t",row.names=F,quote=F)
write.table(counts_table,file=paste0(root_dir,"/stats/SRL_counts_table_n",mismatch,".txt"),sep="\t",row.names=F,quote=F)
