#!/bin/bash

#for loop for defining and finding names of files
for k in {1..10};
do

name=$(sed -n $k{p} ../output/intermediate/51-60.txt)

#creates sam files
bwa mem ../cleaned-data/Septoria/Reference_indexing/ref.fasta ../data/Septoria/reads/*$name*R1*.gz ../data/Septoria/reads/*$name*R2*.gz > ../output/intermediate/$name.sam

#adds read group to sam files
java -jar ../../../programs/picard.jar AddOrReplaceReadGroups -I ../output/intermediate/$name.sam -O ../output/intermediate/$name.rg.sam -RGID $name -RGLB $name -RGPL Illumna -RGPU $name -RGSM $name

#sort the sam file, conver to bam file, and index output
java -jar ../../../programs/picard.jar SortSam -I ../output/intermediate/$name.rg.sam -O ../output/intermediate/$name.bam -SORT_ORDER coordinate -CREATE_INDEX true

#Mark duplicates, write them out to .mdup, filter bam file, and index output
java -jar ../../../programs/picard.jar MarkDuplicates -I ../output/intermediate/$name.bam -O ../output/intermediate/$name.mdup.bam -M ../output/intermediate/$name.mdup -ASSUME_SORT_ORDER coordinate -CREATE_INDEX true

#Sort the marked up bam file and index output
java -jar ../../../programs/picard.jar SortSam -I ../output/intermediate/$name.mdup.bam -O ../output/intermediate/$name.sorted.bam -SORT_ORDER coordinate -CREATE_INDEX true

#produce vcf files from the .sorted.bam files
gatk HaplotypeCaller -R ../cleaned-data/Septoria/Reference_indexing/ref.fasta -I ../output/intermediate/$name.sorted.bam -O ../output/intermediate/$name.vcf

##ONLY FOR DEBUGGING/TESTING, moves to a directory
#mv ../output/intermediate/*$name* ../test/

done

