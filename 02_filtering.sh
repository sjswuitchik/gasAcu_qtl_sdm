# SJSW 2018
# Note: currently quite clunky, would be better optomized with for loops

# Prep libraries for quality filtering

mv Lib01.vcf gasAcu_lib01.vcf
mv Lib02.vcf gasAcu_lib02.vcf
mv Lib03.vcf gasAcu_lib03.vcf 

bgzip gasAcu_lib01.vcf -c > gasAcu_lib01.vcf.gz
tabix -p vcf gasAcu_lib01.vcf.gz 
bgzip gasAcu_lib02.vcf -c > gasAcu_lib02.vcf.gz
tabix -p vcf gasAcu_lib02.vcf.gz
bgzip gasAcu_lib03.vcf -c > gasAcu_lib03.vcf.gz
tabix -p vcf gasAcu_lib03.vcf.gz

bcftools merge gasAcu_lib01.vcf.gz gasAcu_lib02.vcf.gz gasAcu_lib03.vcf.gz -O z -o gasAcu_lib123.vcf.gz
tabix -p vcf gasAcu_lib123.vcf.gz

# Iterative filtering, adapted from dDocent best practices by Jon Puritz
vcftools --gzvcf gasAcu_lib123.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.75 --minGQ 25 --minDP 3 --recode --recode-INFO-all --out gasAcu_lib123_f1
vcftools --vcf gasAcu_lib123_f1.recode.vcf --maf 0.01 --min-meanDP 15 --max-meanDP 200 --recode --recode-INFO-all --out gasAcu_lib123_final
mv gasAcu_lib123_final.recode.vcf gasAcu_lib123_final.vcf
vcftools --vcf gasAcu_lib123_final.vcf --012 --out gasAcu_lib123

# transfer 012.indv list to local machine
# create lists for each of indv to remove from each 'final' vcf on cluster
# took out all Hotel individuals, along with Klein individuals that were not used in this exp

vcftools --vcf gasAcu_lib123_final.vcf --remove hotelindv_toremove.indv --recode --recode-INFO-all --out gasAcu_lib123_klein
mv gasAcu_lib123_klein.recode.vcf gasAcu_lib123_klein.vcf
vcftools --vcf gasAcu_lib123_klein.vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out gasAcu_lib123_oct4
mv gasAcu_lib123_oct4.recode.vcf gasAcu_lib123_oct4.vcf
