gatk --java-options "-Xmx30g -Djava.io.tmpdir=./tmp" MergeVcfs -I raw_vcf.list -O all.merge_raw.vcf
