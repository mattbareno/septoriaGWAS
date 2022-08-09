#!/bin/bash

#remember to enter the cutadapt environment with this command:
        #conda activate cutadaptenv

#run this file in the bash_scripts directory

#creates files of all fastq files in a directory
#ls ../data/*.gz > allfiles.txt
#cuts the file list by unique identifiers
<<<<<<< HEAD
#cut -d "-" -f 7-9 allfiles.txt > dashdelimit.txt
#cut -d "_" -f 1-2 dashdelimit.txt > identifier.txt
=======
ls *.fastq > allfiles.txt
cut -d "-" -f 7-9 allfiles.txt > dashdelimit.txt
cut -d "_" -f 1-2 dashdelimit.txt > identifier.txt
>>>>>>> 3fa1843bdd8d1051af67d8942f6de94229677a35

#for loop for defining and finding names of files
for k in {1..20}; 
do

name=$(sed -n $k{p} identifier.txt)

#cuts primers using paired-end reads, creates intermediate primer-cut fastqs
cutadapt -j 20 -a ATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCTCGTATG -A ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTC -o intermediate$name.R1i.fastq.gz  -p intermediate$name.R2i.fastq.gz ../data/*$name*R1*.gz ../data/*$name*R2*.gz

#filters intermediate primer-cut fastqs by quality >=30. final cleaed output files are in the form UNIQUEIDENTIFIERR(1/2)_001.fastq
cutadapt -j 20 -q 30 -o ../cleaned-data/quality_filtered/$name.R1_001.fastq.gz -p ../cleaned-data/quality_filtered/$name.R2_001.fastq.gz intermediate$name.R1i.fastq.gz intermediate$name.R2i.fastq.gz

#performing fastqc on primer-cut, quality filtered file
fastqc  ../cleaned-data/quality_filtered/$name.R1_001.fastq.gz
fastqc  ../cleaned-data/quality_filtered/$name.R2_001.fastq.gz

done

#moves the fastqc htmls to another folder

mv ../cleaned-data/quality_filtered/*.html ../output/quality_filtered_fastqc_htmls/

#removes intermediate files for storage
rm -r intermediate*.fastq
rm -r ../cleaned-data/quality_filtered/*fastqc.zip
rm -r allfiles.txt
rm -r dashdelimit.txt
rm -r identifier.txt
	
#optional, removes the primer/quality filtered fastq files
#rm -r  ../output/quality_filtered_fastqc_htmls/*R1.fastq
#rm -r  ../output/quality_filtered_fastqc_htmls/*R2.fastq

#how to edit in different context
##change the range in the for loop to the amount of files
## change the range of the cut function in line 9 according to the location of the identifier of the data set
##change the directories depending on the organization of your code
