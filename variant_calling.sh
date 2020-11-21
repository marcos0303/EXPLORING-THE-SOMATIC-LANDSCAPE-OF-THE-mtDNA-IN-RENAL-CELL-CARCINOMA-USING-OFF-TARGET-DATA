# Variant calling using Mutect2 in paired mode (tumor vs normal sample)

for i in MT_*;do tumor=$(echo $i | cut -f3 -d '.'); normal=$(echo $i | cut -f2 -d '.'); java -jar /home/jlanillos/Disco4tb/usr/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar Mutect2 -R /home/jlanillos/Disco4tb/refGen/GATK_g1k_v37/human_g1k_v37.fasta -I $i/$tumor.tumor.bam -tumor $tumor  -I $i/$normal.germline.bam -normal $normal  -L /home/jlanillos/Disco4tb/refGen/BED_files/MTDNA/mtdna.bed -O $i/$tumor.$normal.m2.paired.vcf;   java -jar /home/jlanillos/Disco4tb/usr/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar FilterMutectCalls -V $i/$tumor.$normal.m2.paired.vcf -O $i/$tumor.$normal.filtered.m2.paired.vcf -R /home/jlanillos/Disco4tb/refGen/GATK_g1k_v37/human_g1k_v37.fasta; done
# Compress into -gz all vcf and indexed them (necessary for the multisample merging step)
for i in MT_*; do bgzip $i/*filtered.m2.paired.vcf; tabix -p vcf 4i/*filtered.m2.paired.vcf.gz; done


# Merge all single sample VCF into a multisample VCF

realpath MT_*/*filtered.m2.paired.vcf.gz >> list_singleVCFs.txt
bcftools merge -o multisample.vcf.gz -Oz -l list_singleVCFs.txt

# Variant annotation using the Variant Effect Predictor
conda deactivate; /home/jlanillos/Disco4tb/usr/ensembl-vep/vep -i multisample.vcf.gz -o vep_multisample.vcf --fork 16 --vcf --poly p --sift p --variant_class -format vcf --offline -v --force_overwrite --assembly GRCh38 --everything
conda activate py36

# Split all alleles
bgzip vep_multisample.vcf; tabix -p vcf vep_multisample.vcf.gz
bcftools norm -m -both -o vep_multisample_multiallelic.vcf -Ov VEP_MERGED.vcf.gz

# Use a Python script to parse "vep_multisample_multiallelic.vcf" into a data frame (pending)
