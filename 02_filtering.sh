# SJSW 2018
# Note: currently quite clunky, would be better optomized with for loops

#conda create -n phd -c bioconda plink vcftools htslib bcftools rename
source activate phd 

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

for file in *.gz;
do
  vcftools --gzvcf $file --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.75 --minGQ 25 --minDP 3 --recode --recode-INFO-all --out $file.tmp 2> $file.filter1.log
  vcftools --vcf $file.tmp.recode.vcf --maf 0.01 --min-meanDP 15 --max-meanDP 200 --recode --recode-INFO-all --out $file.filter 2> $file.filter2.log
done

for file in *.filter.recode.vcf;
do
  bgzip $file
  tabix -p vcf $file.gz
done

rm *.tmp.recode.vcf
bcftools merge gasAcu_lib01.vcf.gz.filter.recode.vcf.gz gasAcu_lib02.vcf.gz.filter.recode.vcf.gz gasAcu_lib03.vcf.gz.filter.recode.vcf.gz -O z -o gasAcu_merge.vcf.gz
tabix -p vcf gasAcu_merge.vcf.gz

cp lib123_hotelindv_toremove.indv lib123_kleinindv_toremove.indv /n/holyscratch01/informatics/swuitchik/gasAcu

vcftools --gzvcf gasAcu_merge.vcf.gz --remove lib123_hotelindv_toremove.indv --min-alleles 2 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out gasAcu.merge.kelin 2> final.klein.log

vcftools --gzvcf gasAcu_merge.vcf.gz --remove lib123_kleinindv_toremove.indv --min-alleles 2 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out gasAcu.merge.hotel 2> final.hotel.log
