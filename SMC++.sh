#ln -s  ../data/Chr01.fa  genome.fa

## Generate 35bp sliding-window reads from the genome and split the file
# splitfa genome.fa  35  > read.fa
split -l 20000000 -d   read.fa   read.split

## Build BWA index for the genome
bwa index  genome.fa

## Generate BWA aln alignment scripts
ls ./read.split*[0-9] | while read aa
do
    echo "
bwa aln -t 4  -R 1000000 -O 3 -E 3 genome.fa  $aa > $aa.sai
bwa samse genome.fa  $aa.sai  $aa |gzip > $aa.sam.gz
" > $aa.bwa.sh

done


## Run alignment scripts
# Note: I corrected 'read file' to 'read aa' below to match the '$aa' variable used inside the loop
ls read.split*.bwa.sh | while read aa
do
    sh $aa 1>$aa.log  2>$aa.err
done

## Process alignment results and convert to rawMask file
gzip -dc read.split*.sam.gz | gen_raw_mask.pl > rawMask_35.fa

## Set filtering threshold r=0.5
gen_mask -l 35  -r 0.5 rawMask_35.fa > mask_35_50.fa

## Generate the masked genome file
apply_mask_s mask_35_50.fa  genome.fa  > genome.mask.fa

## Extract masked regions to generate a BED format file
perl  get_lcase.pl  genome.mask.fa  > genome.mask.bed

## Compress and generate index
bgzip  genome.mask.bed
tabix  genome.mask.bed.gz

## Prepare sample information, generating format pop:sp1,sp2
cat ../data/sample.list | awk '$2=="Msie" { sp = sp $1 ","}END{print "Msie:" sp } ' | sed 's/,$//' > Msie.list

## Prepare Composite Likelihood sample list, generally selecting 2-10 samples
cat ../data/sample.list | awk '$2=="Msie"{print $1}' | head -n 4 > Msie_sample.list

vcf=../data/sample.vcf.gz
mask=genome.mask.bed.gz

## Convert VCF to SMC input format
mkdir  Msie_vcf2smc

cat chr.list |while read chr
do
        cat  Msie_sample.list  | while read sample
        do
                # Note: 'smc++' is now called directly without singularity
                echo "smc++ vcf2smc -m $mask -d $sample $sample $vcf  Msie_vcf2smc/$sample.$chr.smc.gz  $chr `cat Msie.list` "
        done
done > Msie_vcf2smc.sh

sh  Msie_vcf2smc.sh

## Estimate historical effective population size
mkdir Msie_analysis

smc++ estimate --cores 10 --knots 10 --spline cubic  -o Msie_analysis  2.9e-8   Msie_v