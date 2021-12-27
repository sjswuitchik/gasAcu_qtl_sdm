# AP 2018

#### STEP1: PROCESS RADTAGS
## 1.1 without any filtering

# Library 1
process_radtags -P -p /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/raw -b /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/raw/barcodes.txt –o /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/raw --renz_1 nlaIII --renz_2 mluCI -i gzfastq -y fastq --inline_index --disable_rad_check -r

# Library 2

process_radtags -P -p /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/raw/raw1 -b /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/raw/raw1/barcodes.txt –o /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/raw/raw1 --renz_1 nlaIII --renz_2 mluCI -i gzfastq -y fastq --inline_index --disable_rad_check -r
process_radtags -P -p /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/raw/raw2 -b /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/raw/raw2/barcodes.txt –o /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/raw/raw2 --renz_1 nlaIII --renz_2 mluCI -i gzfastq -y fastq --inline_index --disable_rad_check -r

# Library 3
process_radtags -P -p /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/raw -b /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/raw/barcodes.txt –o /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/raw --renz_1 nlaIII --renz_2 mluCI -i gzfastq -y fastq --inline_index --disable_rad_check -r

rm *rem.*

## 1.2 Extract fastqc reports

for R1 in *1.fq    
do
name=`echo $R1 | sed 's/.1.fq\+//'`       
R2=${R1%1.fq}2.fq
echo fastqc $R1
echo fastqc $R2
done


# Lib1
mv *1.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/demulti
mv *2.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/demulti
mv *.html /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/demulti/fastqc
mv *.zip /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/demulti/fastqc

multiqc /sf1/project/xpa-194-aa/antoine/temp-tol/Lib01/demulti/fastqc

# Lib2
#1
mv *1.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti1
mv *2.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti1
mv *.html /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti1/fastqc
mv *.zip /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti1/fastqc

multiqc /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti1/fastqc

#2
mv *1.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti2
mv *2.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti2
mv *.html /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti2/fastqc
mv *.zip /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti2/fastqc

multiqc /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/demulti/demulti2/fastqc


# Lib3
mv *1.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/demulti
mv *2.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/demulti
mv *.html /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/demulti/fastqc
mv *.zip /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/demulti/fastqc

multiqc /sf1/project/xpa-194-aa/antoine/temp-tol/Lib03/demulti/fastqc




#### STEP 2: REMOVE PCR DUPLICATES

for R1 in *1.fq
do
name=`echo $R1 | sed 's/.1.fq\+//'`
R2=${R1%1.fq}2.fq       #define a new variable R2 which takes variable R1 and finds (%) 1.fq and replaced with 2.fq
clone_filter -1 $R1 -2 $R2 -i fastq -D --oligo_len_2 7 --index_inline >> $name.log
done



#### STEP 3: HOUSECLEANING I - rename files and concatenate Lib02 files

cp *.1.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/nodup/allnodup/
cp *.2.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/nodup/allnodup/

rename .1.1.fq .1.1.1.fq *.1.1.fq 
rename .2.2.fq .2.2.2.fq *.2.2.fq

cp *.1.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/nodup/allnodup/
cp *.2.fq /sf1/project/xpa-194-aa/antoine/temp-tol/Lib02/nodup/allnodup/

cat RFBxHL17.1.1.fq RFBxHL17.1.1.1.fq > concat_RFBxHL17.1.1.fq
cat RFBxHL17.2.2.fq RFBxHL17.2.2.2.fq > concat_RFBxHL17.1.2.fq
cat RFBxHL18.1.1.fq RFBxHL18.1.1.1.fq > concat_RFBxHL18.1.1.fq
cat RFBxHL18.2.2.fq RFBxHL18.2.2.2.fq > concat_RFBxHL18.1.2.fq
cat RFBxHL19.1.1.fq RFBxHL19.1.1.1.fq > concat_RFBxHL19.1.1.fq
cat RFBxHL19.2.2.fq RFBxHL19.2.2.2.fq > concat_RFBxHL19.1.2.fq
cat RFBxHL20.1.1.fq RFBxHL20.1.1.1.fq > concat_RFBxHL20.1.1.fq
cat RFBxHL20.2.2.fq RFBxHL20.2.2.2.fq > concat_RFBxHL20.1.2.fq
cat RFBxHL21.1.1.fq RFBxHL21.1.1.1.fq > concat_RFBxHL21.1.1.fq
cat RFBxHL21.2.2.fq RFBxHL21.2.2.2.fq > concat_RFBxHL21.1.2.fq
cat RFBxHL22.1.1.fq RFBxHL22.1.1.1.fq > concat_RFBxHL22.1.1.fq 
cat RFBxHL22.2.2.fq RFBxHL22.2.2.2.fq > concat_RFBxHL22.1.2.fq
cat RFBxHL23.1.1.fq RFBxHL23.1.1.1.fq > concat_RFBxHL23.1.1.fq
cat RFBxHL23.2.2.fq RFBxHL23.2.2.2.fq > concat_RFBxHL23.1.2.fq
cat RFBxHL24.1.1.fq RFBxHL24.1.1.1.fq > concat_RFBxHL24.1.1.fq
cat RFBxHL24.2.2.fq RFBxHL24.2.2.2.fq > concat_RFBxHL24.1.2.fq
cat RFBxHL25.1.1.fq RFBxHL25.1.1.1.fq > concat_RFBxHL25.1.1.fq
cat RFBxHL25.2.2.fq RFBxHL25.2.2.2.fq > concat_RFBxHL25.1.2.fq
cat RFBxHL26.1.1.fq RFBxHL26.1.1.1.fq > concat_RFBxHL26.1.1.fq
cat RFBxHL26.2.2.fq RFBxHL26.2.2.2.fq > concat_RFBxHL26.1.2.fq 
cat RFBxHL27.1.1.fq RFBxHL27.1.1.1.fq > concat_RFBxHL27.1.1.fq 
cat RFBxHL27.2.2.fq RFBxHL27.2.2.2.fq > concat_RFBxHL27.1.2.fq 
cat RFBxHL28.1.1.fq RFBxHL28.1.1.1.fq > concat_RFBxHL28.1.1.fq 
cat RFBxHL28.2.2.fq RFBxHL28.2.2.2.fq > concat_RFBxHL28.1.2.fq
cat RFBxHL29.1.1.fq RFBxHL29.1.1.1.fq > concat_RFBxHL29.1.1.fq
cat RFBxHL29.2.2.fq RFBxHL29.2.2.2.fq > concat_RFBxHL29.1.2.fq
cat RFBxHL30.1.1.fq RFBxHL30.1.1.1.fq > concat_RFBxHL30.1.1.fq
cat RFBxHL30.2.2.fq RFBxHL30.2.2.2.fq > concat_RFBxHL30.1.2.fq
cat RFBxKL10.1.1.fq RFBxKL10.1.1.1.fq > concat_RFBxKL10.1.1.fq
cat RFBxKL10.2.2.fq RFBxKL10.2.2.2.fq > concat_RFBxKL10.1.2.fq
cat RFBxKL11.1.1.fq RFBxKL11.1.1.1.fq > concat_RFBxKL11.1.1.fq
cat RFBxKL1.1.1.fq RFBxKL1.1.1.1.fq > concat_RFBxKL1.1.1.fq
cat RFBxKL11.2.2.fq RFBxKL11.2.2.2.fq > concat_RFBxKL11.1.2.fq
cat RFBxKL12.1.1.fq RFBxKL12.1.1.1.fq > concat_RFBxKL12.1.1.fq
cat RFBxKL12.2.2.fq RFBxKL12.2.2.2.fq > concat_RFBxKL12.1.2.fq
cat RFBxKL1.2.2.fq RFBxKL1.2.2.2.fq > concat_RFBxKL1.1.2.fq
cat RFBxKL13.1.1.fq RFBxKL13.1.1.1.fq > concat_RFBxKL13.1.1.fq
cat RFBxKL13.2.2.fq RFBxKL13.2.2.2.fq > concat_RFBxKL13.1.2.fq
cat RFBxKL14.1.1.fq RFBxKL14.1.1.1.fq > concat_RFBxKL14.1.1.fq
cat RFBxKL14.2.2.fq RFBxKL14.2.2.2.fq > concat_RFBxKL14.1.2.fq
cat RFBxKL15.1.1.fq RFBxKL15.1.1.1.fq > concat_RFBxKL15.1.1.fq
cat RFBxKL15.2.2.fq RFBxKL15.2.2.2.fq > concat_RFBxKL15.1.2.fq
cat RFBxKL16.1.1.fq RFBxKL16.1.1.1.fq > concat_RFBxKL16.1.1.fq
cat RFBxKL16.2.2.fq RFBxKL16.2.2.2.fq > concat_RFBxKL16.1.2.fq
cat RFBxKL17.1.1.fq RFBxKL17.1.1.1.fq > concat_RFBxKL17.1.1.fq
cat RFBxKL17.2.2.fq RFBxKL17.2.2.2.fq > concat_RFBxKL17.1.2.fq
cat RFBxKL18.1.1.fq RFBxKL18.1.1.1.fq > concat_RFBxKL18.1.1.fq
cat RFBxKL18.2.2.fq RFBxKL18.2.2.2.fq > concat_RFBxKL18.1.2.fq
cat RFBxKL19.1.1.fq RFBxKL19.1.1.1.fq > concat_RFBxKL19.1.1.fq
cat RFBxKL19.2.2.fq RFBxKL19.2.2.2.fq > concat_RFBxKL19.1.2.fq
cat RFBxKL20.1.1.fq RFBxKL20.1.1.1.fq > concat_RFBxKL20.1.1.fq
cat RFBxKL20.2.2.fq RFBxKL20.2.2.2.fq > concat_RFBxKL20.1.2.fq
cat RFBxKL2.1.1.fq RFBxKL2.1.1.1.fq > concat_RFBxKL2.1.1.fq
cat RFBxKL22.1.1.fq RFBxKL22.1.1.1.fq > concat_RFBxKL22.1.1.fq
cat RFBxKL22.2.2.fq RFBxKL22.2.2.2.fq > concat_RFBxKL22.2.2.fq
cat RFBxKL2.2.2.fq RFBxKL2.2.2.2.fq > concat_RFBxKL2.2.2.fq
cat RFBxKL24.1.1.fq RFBxKL24.1.1.1.fq > concat_RFBxKL24.1.1.fq
cat RFBxKL24.2.2.fq RFBxKL24.2.2.2.fq > concat_RFBxKL24.1.2.fq
cat RFBxKL25.1.1.fq RFBxKL25.1.1.1.fq > concat_RFBxKL25.1.1.fq
cat RFBxKL25.2.2.fq RFBxKL25.2.2.2.fq > concat_RFBxKL25.1.2.fq
cat RFBxKL26.1.1.fq RFBxKL26.1.1.1.fq > concat_RFBxKL26.1.1.fq
cat RFBxKL26.2.2.fq RFBxKL26.2.2.2.fq > concat_RFBxKL26.1.2.fq
cat RFBxKL28.1.1.fq RFBxKL28.1.1.1.fq > concat_RFBxKL28.1.1.fq
cat RFBxKL28.2.2.fq RFBxKL28.2.2.2.fq > concat_RFBxKL28.1.2.fq
cat RFBxKL29.1.1.fq RFBxKL29.1.1.1.fq > concat_RFBxKL29.1.1.fq
cat RFBxKL29.2.2.fq RFBxKL29.2.2.2.fq > concat_RFBxKL29.1.2.fq
cat RFBxKL30.1.1.fq RFBxKL30.1.1.1.fq > concat_RFBxKL30.1.1.fq
cat RFBxKL30.2.2.fq RFBxKL30.2.2.2.fq > concat_RFBxKL30.1.2.fq
cat RFBxKL31.1.1.fq RFBxKL31.1.1.1.fq > concat_RFBxKL31.1.1.fq
cat RFBxKL3.1.1.fq RFBxKL3.1.1.1.fq > concat_RFBxKL3.1.1.fq
cat RFBxKL31.2.2.fq RFBxKL31.2.2.2.fq > concat_RFBxKL31.1.2.fq
cat RFBxKL32.1.1.fq RFBxKL32.1.1.1.fq > concat_RFBxKL32.1.1.fq
cat RFBxKL32.2.2.fq RFBxKL32.2.2.2.fq > concat_RFBxKL32.1.2.fq
cat RFBxKL3.2.2.fq RFBxKL3.2.2.2.fq > concat_RFBxKL3.1.2.fq
cat RFBxKL33.1.1.fq RFBxKL33.1.1.1.fq > concat_RFBxKL33.1.1.fq
cat RFBxKL33.2.2.fq RFBxKL33.2.2.2.fq > concat_RFBxKL33.1.2.fq
cat RFBxKL34.1.1.fq RFBxKL34.1.1.1.fq > concat_RFBxKL34.1.1.fq
cat RFBxKL34.2.2.fq RFBxKL34.2.2.2.fq > concat_RFBxKL34.1.2.fq
cat RFBxKL35.1.1.fq RFBxKL35.1.1.1.fq > concat_RFBxKL35.1.1.fq
cat RFBxKL35.2.2.fq RFBxKL35.2.2.2.fq > concat_RFBxKL35.1.2.fq
cat RFBxKL36.1.1.fq RFBxKL36.1.1.1.fq > concat_RFBxKL36.1.1.fq
cat RFBxKL36.2.2.fq RFBxKL36.2.2.2.fq > concat_RFBxKL36.1.2.fq
cat RFBxKL38.1.1.fq RFBxKL38.1.1.1.fq > concat_RFBxKL38.1.1.fq
cat RFBxKL38.2.2.fq RFBxKL38.2.2.2.fq > concat_RFBxKL38.1.2.fq
cat RFBxKL4.1.1.fq RFBxKL4.1.1.1.fq > concat_RFBxKL4.1.1.fq
cat RFBxKL4.2.2.fq RFBxKL4.2.2.2.fq > concat_RFBxKL4.1.2.fq
cat RFBxKL5.1.1.fq RFBxKL5.1.1.1.fq > concat_RFBxKL5.1.1.fq
cat RFBxKL5.2.2.fq RFBxKL5.2.2.2.fq > concat_RFBxKL5.1.2.fq
cat RFBxKL6.1.1.fq RFBxKL6.1.1.1.fq > concat_RFBxKL6.1.1.fq
cat RFBxKL6.2.2.fq RFBxKL6.2.2.2.fq > concat_RFBxKL6.1.2.fq
cat RFBxKL7.1.1.fq RFBxKL7.1.1.1.fq > concat_RFBxKL7.1.1.fq
cat RFBxKL7.2.2.fq RFBxKL7.2.2.2.fq > concat_RFBxKL7.1.2.fq
cat RFBxKL8.1.1.fq RFBxKL8.1.1.1.fq > concat_RFBxKL8.1.1.fq
cat RFBxKL8.2.2.fq RFBxKL8.2.2.2.fq > concat_RFBxKL8.1.2.fq
cat RFBxKL9.1.1.fq RFBxKL9.1.1.1.fq > concat_RFBxKL9.1.1.fq
cat RFBxKL9.2.2.fq RFBxKL9.2.2.2.fq > concat_RFBxKL9.1.2.fq




#### STEP4: ALIGNMENT BWA mem

rename .2.2.fq .1.2.fq *.2.2.fq


for R1 in *1.fq    
do
name=`echo $R1 | sed 's/.1.fq\+//'`       
R2=${R1%1.fq}2.fq
echo "bwa mem -t 16 -R '@RG\tID:"$name'\tPL:Illumina\tLB:'$name'\tSM:'$name"'" /sf1/project/xpa-194-aa/stickl_genome/stickl_genome85/stickl_BWA/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa $R1 $R2 '>' $name.sam '2>' $name.BWAmem.log
done




#### STEP5: Remove ambiguous reads from sam files and make sorted.bam files


for i in *.sam  
 do  
   echo $i  
   name=$(echo ${i} | sed 's/.sam//')
   echo $name  

   samtools view -q 20 -b -S $i > ${name}.bam
   samtools sort ${name}.bam -o ${name}_sorted.bam

 done




#### STEP6: CREATE A BAM LIST

for f in `ls -S *_sorted.bam`
do echo $f >> bam_list.txt
done

chmod a+rx bam_list.txt



#### STEP7: SAMTOOLS MPILEUP

samtools mpileup -C 50 -E -t SP -t DP -u -I -f /sf1/project/xpa-194-aa/stickl_genome/stickl_genome85/stickl_samtools/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa -b bam_list.txt > Lib01.bcf
samtools mpileup -C 50 -E -t SP -t DP -u -I -f /sf1/project/xpa-194-aa/stickl_genome/stickl_genome85/stickl_samtools/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa -b bam_list.txt > Lib02.bcf
samtools mpileup -C 50 -E -t SP -t DP -u -I -f /sf1/project/xpa-194-aa/stickl_genome/stickl_genome85/stickl_samtools/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa -b bam_list.txt > Lib03.bcf



#### STEP8: MAKE VCF FILE

bcftools call -v -c -f GQ Lib01.bcf > Lib01.vcf
bcftools call -v -c -f GQ Lib02.bcf > Lib02.vcf
bcftools call -v -c -f GQ Lib03.bcf > Lib03.vcf
