#!/usr/bin/bash

name=$1

gatk HaplotypeCaller -R ../cleaned-data/Septoria/Reference_indexing/ref.fasta -I ../output/intermediate/$name.sorted.bam -O ../output/intermediate2/$name.vcf -ERC GVCF -ploidy 1
