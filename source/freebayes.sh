#! /bin/bash
# $1=sample

 . /etc/profile.d/modules.sh

module load briglo/freebayes/v1.2.0-4-gd15209e  briglo/sc_split/190313 pethum/vcftools/gcc-4.4.6/0.1.15

hdir="/share/ScratchGeneral/briglo/scRNA/venchi"
mkdir "$hdir"/outputs/freebayes/"$1"/
###########DONE##############
#for i in * ; do qsub -pe smp 4 -b y -j y -N rmd -V -cwd "samtools view -h "$i"/outs/possorted_genome_bam.bam hg19_1 hg19_2 hg19_3 hg19_4 hg19_5 hg19_6 hg19_7 hg19_8 hg19_9 hg19_10 hg19_11 hg19_12 hg19_13 hg19_14 hg19_15 hg19_16 hg19_17 hg19_18 hg19_19 hg19_20 hg19_21 hg19_22 | samtools view -h -S -q 10 -F 3844 - | samtools sort -@ 4 -m 4G -o "$i"/outs/human_trim_sort_"$i".bam -" ; done
#for i in * ; do qsub -b y -j y -N rmd -V -cwd samtools rmdup -s --reference /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa "$i"/outs/human_trim_sort_"$i".bam "$i"/outs/rmdup_human_trim_sort_"$i".bam ; done 
#for i in * ; do qsub b y -j y -N rmd -V -cwd samtools index "$i"/outs/rmdup_human_trim_sort_"$i".bam ; done 
###extracting "pure" barcode lists for bam subsetting
# md<-read.csv("extraData/pilot_1_metadata.csv",header=T,stringsAsFactors=F)
# x<-split(md,md$orig.ident)
# extBC<-function(metaData){
# bcs<- gsub("_[0-9]",'-1',metaData$Barcode)
# return(list('cancer'=bcs[grep("Cancer",metaData$final.ident)],
#         'normal'=bcs[metaData$final.ident %in% c("Monocyte","CD14+,Macrophage","CD4,T-cell,2","T-Reg","CD8,T-cell,2","Dendritic,cells","NK","CD16+,Macrophage","pDC","CD8,T-cell,1","B-cell","CD4,T-cell,1","Plasmablast")]
#         ))
# }
# bclist<-lapply(x,extBC)
# write.table(bclist$Sample_A$cancer,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleACITE/outs/cancer_barcodes.txt")
# write.table(bclist$Sample_A$normal,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleACITE/outs/normal_barcodes.txt")
# write.table(bclist$Sample_B$cancer,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleBCITE/outs/cancer_barcodes.txt")
# write.table(bclist$Sample_B$normal,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleBCITE/outs/normal_barcodes.txt")

# ### extracting bam records
# module load briglo/subset-bam/1.0_precompiled
# cd /share/ScratchGeneral/briglo/scRNA/venchi/data/hg19/VENCHI_SampleACITE/outs
# qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleACITE.bam --bam-tag CB:Z --cell-barcodes normal_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleACITE_normal.bam
# qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleACITE.bam --bam-tag CB:Z --cell-barcodes cancer_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleACITE_cancer.bam
# samtools index rmdup_human_trim_sort_VENCHI_SampleACITE_normal.bam
# samtools index rmdup_human_trim_sort_VENCHI_SampleACITE_cancer.bam
# cd /share/ScratchGeneral/briglo/scRNA/venchi/data/hg19/VENCHI_SampleBCITE/outs
# qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleBCITE.bam --bam-tag CB:Z --cell-barcodes normal_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleBCITE_normal.bam
# qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleBCITE.bam --bam-tag CB:Z --cell-barcodes cancer_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleBCITE_cancer.bam
# samtools index rmdup_human_trim_sort_VENCHI_SampleBCITE_normal.bam
# samtools index rmdup_human_trim_sort_VENCHI_SampleBCITE_cancer.bam

### new strategy for freebayes 31/5/19 rename readgroups in bams and combine need to resort to scatter gather freebayes
module load briglo/freebayes/v1.2.0-4-gd15209e 
cd /share/ScratchGeneral/briglo/scRNA/venchi/data/hg19/VENCHI_SampleACITE
 qsub -V -cwd -pe smp 2 -N renMerge_A -m e -M b.gloss@garvan.org.au -b y -j y "bamaddrg -b rmdup_human_trim_sort_VENCHI_SampleACITE_cancer.bam -s ACITE_cancer -b rmdup_human_trim_sort_VENCHI_SampleACITE_normal.bam -s ACITE_normal | samtools sort -o ACITE_MergeRG.bam -"
   cd ../../VENCHI_SampleBCITE/outs
   qsub -V -cwd -pe smp 2 -N renMerge_B -m e -M b.gloss@garvan.org.au -b y -j y "bamaddrg -b rmdup_human_trim_sort_VENCHI_SampleBCITE_cancer.bam -s BCITE_cancer -b rmdup_human_trim_sort_VENCHI_SampleBCITE_normal.bam -s BCITE_normal | samtools sort -o BCITE_MergeRG.bam -"


###running freebayes
hdir="/share/ScratchGeneral/briglo/scRNA/venchi"
module load briglo/freebayes/v1.2.0-4-gd15209e 
for i in  `ls data/hg19/` ; do mkdir "$hdir"/outputs/freebayes/"$i"/ 
for j in `cat "$hdir"/source/chrom.txt` ; do \
qsub -pe smp 2 -b y -j y -N fb_"$i"_"$j" -V -cwd "samtools view -h "$hdir"/data/hg19/"$i"/outs/rmdup_human_trim_sort_"$i"_normal.bam "$j" | freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 --min-coverage 12 --stdin  -v "$hdir"/outputs/freebayes/"$i"/"$i"_"$j"_normal_SNV.vcf"
qsub -pe smp 2 -b y -j y -N fb_"$i"_"$j" -V -cwd "samtools view -h "$hdir"/data/hg19/"$i"/outs/rmdup_human_trim_sort_"$i"_cancer.bam "$j" | freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 --min-coverage 12 --stdin -v "$hdir"/outputs/freebayes/"$i"/"$i"_"$j"_cancer_SNV.vcf"
done
done
###compress
module load briglo/samtools/1.5
for i in *vcf ; do qsub -V -cwd -b y -j y bgzip -i "$i" ; done
###combine
module load pethum/vcftools/gcc-4.4.6/0.1.15
ls *normal*.gz > norm.txt
ls *cancer*.gz > can.txt
vcf-concat -f norm.txt > SampleBCITE_normal.vcf.gz
vcf-concat -f can.txt > SampleBCITE_cancer.vcf.gz
####QC
module load briglo/bcftools/1.9
bcftools filter -i 'QUAL>30' -O z -o SampleBCITE_cancer_filter.vcf.gz SampleBCITE_cancer.vcf.gz
bcftools filter -i 'QUAL>30' -O z -o SampleBCITE_normal_filter.vcf.gz SampleBCITE_normal.vcf.gz
#############################
 echo "SampleBCITE_cancer" > can
 bcftools reheader --samples can -o SampleBCITE_cancer_filter_rename.vcf.gz SampleBCITE_cancer_filter.vcf.gz
 tabix -p vcf SampleBCITE_normal_filter.vcf.gz
  tabix -p vcf SampleBCITE_cancer_filter_reheader.vcf.gz
  bcftools merge -O z -o SampleBCITE_merged.vcf.gz SampleBCITE_normal_filter.vcf.gz SampleBCITE_cancer_filter_rename.vcf.gz

 echo "SampleACITE_cancer" > can
 bcftools reheader --samples can -o SampleACITE_cancer_filter_rename.vcf.gz SampleACITE_cancer_filter.vcf.gz
 tabix -p vcf SampleACITE_normal_filter.vcf.gz
  tabix -p vcf SampleACITE_cancer_filter_rename.vcf.gz
  
  
  bcftools merge -O z SampleACITE_normal_filter.vcf.gz SampleACITE_cancer_filter_rename.vcf.gz | bcftools sort -O z -o SampleACITE_merged.vcf.gz

  ### rn demuxlet, had to use GT as opposed to GP (does imputation calc this?)
  module load briglo/demuxlet/1.9
  qsub -V -cwd -pe smp 2 -b y -j y -N dm demuxlet --sam ./data/hg19/VENCHI_SampleACITE/outs/rmdup_human_trim_sort_VENCHI_SampleACITE.bam --vcf ./outputs/freebayes/VENCHI_SampleACITE/SampleACITE_merged.vcf.gz --out ./outputs/demuxlet/SampleACITE  --field GT --geno-error 0.001
  qsub -V -cwd -pe smp 2 -b y -j y -N dm demuxlet --sam ./data/hg19/VENCHI_SampleBCITE/outs/rmdup_human_trim_sort_VENCHI_SampleBCITE.bam --vcf ./outputs/freebayes/VENCHI_SampleBCITE/SampleBCITE_merged.vcf.gz --out ./outputs/demuxlet/SampleBCITE  --field GT --geno-error 0.001

###looking in R
fnam<-dir(pattern='best')
  dm<-lapply(fnam, function(x) data.frame(read.table(x,header=T,stringsAsFactors=F),row.names=1)[-1,]) 
  md<-read.csv("../../extraData/pilot_1_metadata.csv",header=T,stringsAsFactors=F)
rownames(md)<-md$Barcode
