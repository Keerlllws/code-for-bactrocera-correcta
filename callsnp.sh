gatk  --java-options "-Xmx4g -Djava.io.tmpdir=./tmp"  SelectVariants  -R /home/lws2/workspace/BDYT/genome.fasta -V Chr01.raw.vcf.gz --select-type SNP -O Chr01.raw.snp.vcf

gatk  --java-options "-Xmx4g -Djava.io.tmpdir=./tmp"  VariantFiltration -R /home/lws2/workspace/BDYT/genome.fasta -V Chr01.raw.snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter' -O Chr01.filter.snp.vcf

gatk  --java-options "-Xmx4g -Djava.io.tmpdir=./tmp"  SelectVariants  -R /home/lws2/workspace/BDYT/genome.fasta -V Chr01.filter.snp.vcf --exclude-filtered  -O Chr01.filtered.snp.vcf
