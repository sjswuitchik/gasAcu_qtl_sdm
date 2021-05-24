# SJSW 2021
# Note: currently quite clunky and repetitive, would be better optomized with for loops

source activate phd

# Download snpEff, which includes a useful splitting function to separate a single vcf into chromosomal vcfs

mkdir chrSplit
mkdir -p lepmap/hotel
cd chrSplit
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
mkdir -p snpEff/data/hotel
cp ../gasAcu.merge.klein.recode.vcf snpEff/data
cp ../gasAcu.merge.hotel.recode.vcf snpEff/data/hotel
cd snpEff
java -jar SnpSift.jar split data/gasAcu.merge.klein.recode.vcf
java -jar SnpSift.jar split data/hotel/gasAcu.merge.hotel.recode.vcf
cp -v data/gasAcu.merge.klein.recode.group* ../../lepmap
cp -v data/hotel/gasAcu.merge.hotel.recode.group* ../../lepmap/hotel
cd ../../lepmap

# Download Lep-MAP3 to local machine then sftp to cluster
wget https://sourceforge.net/projects/lep-map3/files/latest/download
unzip binary+code.zip
cp /n/boslfs02/LABS/informatics/swuitchik/gasAcu/gasAcu_qtl/lepmap/klein/gasAcu_lib123_ped.txt /n/boslfs02/LABS/informatics/swuitchik/gasAcu/gasAcu_qtl/lepmap/gasAcu_lib123_hotel_ped.txt .
mv gasAcu_lib123_hotel_ped.txt hotel/
conda deactivate 

# initialize LM3
module load jdk/16-fasrc01
javac src/*.java -d bin/
java -cp bin/ ParentCall2 data=gasAcu_lib123_ped.txt vcfFile=gasAcu.merge.klein.recode.groupI.vcf halfSibs=1 removeNonInformative=1 > may2021_chr1.call
java -cp bin/ Filtering2 data=may2021_chr1.call removeNonInformative=1 dataTolerance=0.001 > chr1_filter.call
java -cp bin/ SeparateChromosomes2 data=chr1_filter.call lodLimit=6 sizeLimit=5 > chr1_map1.txt
java -cp bin/ JoinSingles2All data=chr1_filter.call map=chr1_map1.txt lodLimit=3 iterate=1 > chr1_map1js.txt
sort chr1_map1js.txt | uniq -c | sort -n
java -cp bin/ OrderMarkers2 data=chr1_filter.call map=chr1_map1js.txt useKosambi=1 minError=0.001 sexAveraged=1 outputPhasedData=1 identicalLimit=0.01 > chr1_ordered.txt

### do for all chrs 
wget https://avikarn.com/image/lepmap/map2genotypes.awk
chmod +x map2genotypes.awk
cp map2genotypes.awk hotel/

for file in *_ordered.txt;
do
	awk -vFullData=1 -f map2genotypes.awk $file > $file.geno
done


# transfer *.geno and may2021_chr1.call (or any of the chr call files) files to local machine, manipulate in Excel to prep for R/qtl load

# Manual prep steps

## create your header - this will be put into every chr geno matrix
# copy row 3 of may2021_chr1.call and insert above geno contents in Excel, starting at the second column - this is your indv header
# label first column 'marker'
# remove all WT and F1 indvs from the header row

## prep a chromosome geno matrix
# open chr ordered geno file, copy everything, and paste into Excel
# delete fourth column (duplicated sex-averaged position)
# replace all values in second column with chr number (e.g. 1)
# insert header above matrix
# find & replace all values in the matrix as: 
#			1 1 = AA
#			1 2 = AB
#			2 1 = BA
#			2 2 = BB
#
# Alternatively, you can replace 11 with A and 22 with B, leaving 12 and 21 to be replaced with H - depends on the crossing scheme used to produce the F2s. It's easiest if you select the 'match case' and 'find entire cells only' options in the find & replace menu to make sure you're not getting spurious replacements in other labels (e.g. sample names)
#
# Then, copy the entire matrix, add another sheet, right click on the top cell of the sheet, paste special -> transpose
#
# remove 'marker', 'CHR', and 'POS' in column 1
# insert as many columns as you need for phenotypes between the first field (indv) and the first marker/chr/pos field
# delete sheet 1 (so it doesn't get read by R when you load the data)
# save as CSV for easier loading with read.cross

# split full geno matrices for all families into family-specific files



### repeat for Hotel from /n/holyscratch01/informatics/swuitchik/gasAcu/lepmap/hotel

