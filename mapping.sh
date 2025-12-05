#!/bin/bash 

ref=/home/data/lws2/pop/genome/genome.fasta

cat sampleCM.txt |while read  sp fq1 fq2 
do 
    echo   " 
bwa mem -t 20 -R '@RG\tID:$sp\tSM:$sp\tPL:illumina'  $ref  $fq1  $fq2  2>$sp.bwa.log | samtools sort -@ 20 -m 10G -o $sp.sort.bam -

/usr/lib/jvm/java-17-openjdk-amd64/bin/java -Xmx8g -XX:ParallelGCThreads=2 -jar /home/data/lws2/pop/software/picard.jar   MarkDuplicates I=$sp.sort.bam O=$sp.sort.markdup.bam  CREATE_INDEX=true  REMOVE_DUPLICATES=true M=$sp.marked_dup_metrics.txt

/pub/miniconda3/bin/samtools  flagstat  $sp.sort.bam > $sp.sort.bam.flagstat

/pub/miniconda3/bin/samtools coverage  $sp.sort.bam > $sp.sort.bam.coverage
" >   mapping.$sp.sh 
done
