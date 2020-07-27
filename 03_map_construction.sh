# SJSW 2019
# Note: currently quite clunky, would be better optomized with for loops

# Download snpEff, which includes a useful splitting function to separate a single vcf into chromosomal vcfs

mkdir chr_split
cd chr_split
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
cd snpEff
java -jar SnpSift.jar split data/gasAcu_lib123_oct4.vcf 
mv gasAcu_lib123_oct4.group* /scratch/swuitchik/gasAcu_qtl/lepmap/klein

# Download Lep-MAP3
wget https://sourceforge.net/projects/lep-map3/files/latest/download
unzip binary+code.zip
# initialize LM3
javac src/*.java -d bin/
# chr1
java -cp bin/ ParentCall2 data=gasAcu_lib123_ped.txt vcfFile=gasAcu_lib123_oct4.groupI.vcf halfSibs=1 removeNonInformative=1 > nov_chr1.call
java -cp bin/ Filtering2 removeNonInformative=1 dataTolerance=0.01 data=nov_chr1.call > chr1_filter.call
java -cp bin/ SeparateChromosomes2 data=chr1_filter.call lodLimit=3 sizeLimit=5 > chr1_map1.txt
java -cp bin/ JoinSingles2All data=chr1_filter.call map=chr1_map1.txt lodLimit=3 iterate=1 > chr1_map1js.txt
sort chr1_map1js.txt | uniq -c | sort -n
java -cp bin/ OrderMarkers2 data=chr1_filter.call map=chr1_map1js.txt useKosambi=1 minError=0.001 sexAveraged=1 outputPhasedData=0 identicalLimit=0.01 > chr1_ordered.txt

## repeat for each chromosome

## cut each map for lm2rqtl_chr.py script
cut -f 1-2 chr1_ordered.txt > marker_pos_chr1.txt

# run each map through lm2rqtl.py script
python lm2rqtl_chr1.py marker_pos_chr1.txt gasAcu_lib123_oct4.groupI.vcf 

plink --vcf gasAcu_lib123_oct4_concatID_chr1.vcf --extract LepMap3_snps_for_vcftools_chr1.txt --allow-extra-chr --recode vcf --out klein_rqtl_chr1
vcftools --vcf klein_rqtl_chr1.vcf --012 --out rqtl_chr1
