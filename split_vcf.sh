bgzip -c myvcf.vcf > myvcf.vcf.gz

tabix -p vcf myvcf.vcf.gz

tabix myvcf.vcf.gz chr1 > chr1.vcf
