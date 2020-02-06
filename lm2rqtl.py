#!/usr/bin/env python/3.6.2
#Lepmap2rqtl.py

#this script takes LepMap3 output files (.txt) and creates two files, one to use for vcftools to keep only the correct snps, and one to use for R/qtl as a header (containing LG, CM position, and marker name).
#Inputs are: 
#(1) output file from LepMap containing all of your linkage groups in one file 
#that has been cut for only the important lines like so:
#   cut -f 1-2 map.txt > marker_pos_chr1.txt
#(2) the original vcf that contained this data
#sample command:
#   python lm2rtqtl_chr1.py marker_pos_chr1.txt gasAcu_lib123_oct4.groupI.vcf

#--------------------------------------------------------##-------------------------------------------------------#
import sys
import re
import pandas as pd
import numpy as np
import os
import io
#--------------------------------------------------------##-------------------------------------------------------#
##FUNCTIONS##

#rewrite last line to contain say "#***", this is so that we can get the last LG as well
def MakeFileReady(infile, outfile):
    lines = open(infile).read().splitlines()
    lines[-1] = '#***'
    open(outfile,'w').write('\n'.join(lines))

#define the function that will be used to extract the LG's for the for loop
def GetLGs(infile):
    with open(infile) as fp:
        LG = re.findall(r'\#.+ (LG = \d+)',fp.read())
    return LG

#define the function that will get the LG information from the file. 
#This function opens the file (outfile from MakeFileReady), looks for the LG that you give it (delim1)
#looks for the end of the LG (represented by '#***'), and returns everything in between
#RE explanation: '\\b{}\\b' looks for exact matches, '(.*?)' captures the first instance of the match
# '[#*][#*]' finds the #*** delimiter, can't use a + because the + matches one or more times
def GetLGfromfile(infile, delim1):
    with open(infile) as fp:
        for LGinfo in re.findall('\\b{}\\b(.*?)[(#*)][(#*)]'.format(delim1), fp.read(), re.S):
            return LGinfo

#to use this function, you will need a data frame named "all_LG" for the split data to be appended to
def SeparateAndCatLGinfo(LGinfo, LG): #lginfo because we're going to split it, LG because we're going to add that as a column
    positions = []
    sites = []
    for line in re.findall('.+\n', LGinfo):
        if re.match("\#", line) or re.match('.+likelihood.+', line): #skips lines that start with either a # or contain likelihood
            pass
        else:
            for result in re.findall('(.+)\t(.+)\n', line): #finds the two values and stores them as a tuple
                sites.append(result[0])
                positions.append(result[1])
    tempframe = pd.DataFrame({'LG': np.repeat(LG+1, len(positions)), 'POS': positions, 'ID': sites}) #makes a dataframe by making a dictionary of your lists
    return tempframe #returns that tempframe to the workspace


def read_vcf_data(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_table(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str,
                'ALT': str,'QUAL': str,'FILTER': str, 'INFO': str}).rename(columns={'#CHROM': 'CHROM'})

def read_vcf_header(path):
    with open(path, 'r') as f:
        header = [l for l in f if l.startswith('##')]
    return header

#--------------------------------------------------------##-------------------------------------------------------#
##FILE HANDLING##

infile = str('marker_pos_chr1.txt') #have to take second argv because the first one is the name of the program
infile_ready = infile.split('.')[0] + '_ready.' + infile.split('.')[1]
infile_vcf = str('gasAcu_lib123_oct4.groupI.vcf')

MakeFileReady(infile, infile_ready)
LG = GetLGs(infile_ready)
all_LG = pd.DataFrame()
for i in range(len(LG)):
    LGinfo = GetLGfromfile(infile_ready, LG[i])
    tempframe = SeparateAndCatLGinfo(LGinfo,i)
    all_LG = all_LG.append(tempframe)
print('GetLGfromfile has found the following linkage groups in this dataset: \n {}'.format([i for i in LG]))
#you now have a dataframe with all of the LG's, Locations, and ID's in it

#Next, let's read in the vcf information
vcf_df = read_vcf_data(infile_vcf)
vcf_header = read_vcf_header(infile_vcf)
#and populate the ID column before renaming the header back to #CHROM
vcf_df['ID'] = vcf_df['CHROM'] + ':' + vcf_df['POS'].map(str)
vcf_df = vcf_df.rename(columns ={'CHROM':'#CHROM'})
#you now have two dataframes that contain information from your original vcf
#and from your LepMap3 output file.
#We'll use these dataframes to write two files.
#one for the r/qtl header
all_LG.to_csv('./LepMap3_header_for_rqtl_chr1.txt', index=False, header = True)
sitenums = [int(i) for i in all_LG['ID']]
sitenames = [str(i) for i in vcf_df['ID']]
sitenums.sort()
#and one for the snps for vcftools
with open('LepMap3_snps_for_vcftools_chr1.txt', 'w') as file:
    for i in range(len(sitenums)):
       file.write(str(sitenames[sitenums[i]])+'\n')

##Finally, we need to make the new VCF that includes the fake ID's:
with open('./{}_concatID_chr1.vcf'.format('gasAcu_lib123_oct4.groupI.vcf'.split('.')[0]), 'w') as file:#sys argv[2]
    for i in range(len(vcf_header)):
        file.write(vcf_header[i])
    vcf_df.to_csv(file, index=False, header = True, sep = '\t')

#--------------------------------------------------------##-------------------------------------------------------#
##CLOSING STATEMENTS##

#remove any temporary files
os.remove(infile_ready)
print("Temporary File(s) Removed!")

#make a list of the files that we've made for the user.
printstring = ("Lepmap2rqtl has created the following files in your current directory: \n" +
               "1) An alternative version of your input vcf " +
               "called {}_concatID_chr1.vcf ".format('gasAcu_lib123_oct4.groupI.vcf'.split('.')[0]) +
               "which has populated the ID column with 'scaffold:position' \n" +
               "2) A file called 'LepMap3_header_for_rqtl_chr1.txt' which can be used" +
               "as your r/qtl dataframe header or to plot your Linkage Groups using r/LinkageMapView \n" +
               "3) A file called 'LepMap3_snps_for_vcftools_chr1.txt' which contains the named snps " +
               "to be used with {}_concatID_chr1.vcf".format('gasAcu_lib123_oct4.groupI.vcf'.split('.')[0]) +
               "to keep only informative sites") 
print(printstring)

#lepmap2rqtl.py