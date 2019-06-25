#! /bin/bash
# $1=sample

 . /etc/profile.d/modules.sh

module load briglo/freebayes/v1.2.0-4-gd15209e  briglo/sc_split/190313 pethum/vcftools/gcc-4.4.6/0.1.15

ddir="/share/ScratchGeneral/briglo/scRNA/chrcha/data"
cd "$ddir"/"$1"/outs

#samtools view -S -b -q 10 -F 3844 "$ddir"/"$1"/outs/possorted_genome_bam.bam | samtools sort -@ 4 -m 4G - | samtools rmdup --reference /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa - "$ddir"/"$1"/outs/rmdup_sort_"$1".bam

  #samtools rmdup "$ddir"/"$1"/outs/"$1"_target.bam "$ddir"/"$1"/outs/rmdups_"$1"_target.bam

 #samtools sort -@ 4 -m 4G -o "$ddir"/"$1"/outs/sort_rmdups_"$1"_target "$ddir"/"$1"/outs/rmdups_"$1"_target.bam #crapping out here, throttled thoughput with other jobs, might actually work

 #samtools index "$ddir"/"$1"/outs/rmdup_sort_"$1".bam

#freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 -v "$ddir"/"$1"/outs/"$1"_SNV.vcf "$ddir"/"$1"/outs/rmdup_sort_"$1".bam

samtools view "$ddir"/"$1"/outs/rmdup_sort_"$1".bam hg19_1 hg19_2 hg19_3 hg19_4 hg19_5 hg19_6 hg19_7 hg19_8 hg19_9 hg19_10 hg19_11 hg19_12 hg19_13 hg19_14 hg19_15 hg19_16 hg19_17 hg19_18 hg19_19 hg19_20 hg19_21 hg19_22 | freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 --min-coverage 12 --stdin -v "$ddir"/"$1"/outs/"$1"_SNV.vcf  -

vcftools --vcf "$ddir"/"$1"/outs/"$1"_SNV.vcf --minQ 30 --recode --recode-INFO-all --out /share/ScratchGeneral/briglo/scRNA/chrcha/vcfs/"$1"_hiqualSNV.vcf"


samtools view --output-fmt bam rmdup_sort_M43973_PT.bam hg19_1 hg19_2 hg19_3 hg19_4 hg19_5 hg19_6 hg19_7 hg19_8 hg19_9 hg19_10 hg19_11 hg19_12 hg19_13 hg19_14 hg19_15 hg19_16 hg19_17 hg19_18 hg19_19 hg19_20 hg19_21 hg19_22 | freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 -v pipe_SNV.vcf -

##### had to make new fa
# sed 's/^>/>hg19_/' /share/ClusterShare/biodata/contrib/briglo/hg19/hs37d5.fa > /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa 
#samtools faidx renamed_hg19.fa
#extract chrom interested in edit in text editor to get lists, parsed to initial samtool view to save on filesize
#grep -e "^>" ../../../anno/renamed_hg19.fa |  cut -d " " -f 1 > chrnames.txt
# liver uses pipe 
#   cd "$ddir"/"$1"/outs
#   samtools view -S -b -q 10 -F 3844 possorted_genome_bam.bam | samtools sort -m 16G - | samtools rmdup - sort_rmdups_"$1"_target (dunno a bit shite)



qsub -V -cwd -l h_vmem=8G -pe smp 2  -N FB_pipe_"$i" -b y -j y "samtools view -O bam "$ddir"/"$i"/outs/rmdup_sort_"$i".bam hg19_1 hg19_2 hg19_3 hg19_4 hg19_5 hg19_6 hg19_7 hg19_8 hg19_9 hg19_10 hg19_11 hg19_12 hg19_13 hg19_14 hg19_15 hg19_16 hg19_17 hg19_18 hg19_19 hg19_20 hg19_21 hg19_22 | freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 --min-coverage 12 --stdin -v "$ddir"/"$i"/outs/"$i"_SNV.vcf"