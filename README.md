# septoria musiva GWAS
STEPS

*this assumes you are given the fastq reads and the reference genome*

1. index the reference fasta

`gatk CreateSequenceDictionary -R REF.fasta`

`samtools faidx REF.fasta`

`bwa index REF.fasta`

*keep everything here in the Reference_indexing directory*

2. combine R1 & R2 of the reads into a SAM file

`bwa mem REF.fasta ID_R1.fastq ID_R2.fastq > ID.sam`

*produces .sam file*

3. add read group to the sam file

`java -jar DIRECTORY/picard.jar AddOrReplaceReadGroups -I ID.sam -O ID.rg.sam -RGID ID -RGLB ID -RGPL Illumna -RGPU ID -RGSM ID`

*produces rg.sam files*

4. sort the sam file, conver to bam file, and index output

`java -jar DIRECTORY/picard.jar SortSam -I ID.rg.sam -O ID.bam -SORT_ORDER coordinate -CREATE_INDEX true`

*produces .bam, .bai files*

5. Mark duplicates, write them out to .mdup, filter bam file, and index output

`java -jar DIRECTORY/picard.jar MarkDuplicates -I ID.bam -O ID.mdup.bam -M ID.mdup -ASSUME_SORT_ORDER coordinate -CREATE_INDEX true`

*produces .mdup.bam, .mdup.bai files*

6. Sort the marked up bam file and index output

`java -jar picard.jar SortSam -I ID.mdup.bam -O ID.sorted.bam -SORT_ORDER coordinate -CREATE_INDEX true`

*produces .sorted.bam, .sorted.bai files*

7. produce vcf files from the .sorted.bam files

`gatk HaplotypeCaller -R REF.fasta -I ID.sorted.bam -O ID.vcf -ERC GVCF -ploidy 1`

*produces .vcf files*

8. produce a combined vcf file using CombineGVCFs


`gatk CombineGVCFs -R REF.fasta -V ID.vcf -V ID2.vcf (...) -V IDN.vcf -O COMBINED.vcf`

a. if you have a lot of files, create a script that prints all of the file names delimited by " -V " 

*produced .vcf FILE*

9. produce a joint genotyping

`gatk GenotypeGVCFs -R REF.fasta -V COMBINED.vcf -O COMBINED.gvcf.vcf`

*produces .vcf FILE*


10. Generate .bed, .bim, .fam files

`plink --vcf X.vcf --allow-extra-chr`

11. Generate .map, .ped, .nosex, and .log files

`plink --vcf X.vcf --recode --out PREFIX --allow-extra-ch`

12. be sure to rename the .bed, .bim. and .fam files to contain the .ped file's name. i.e., if your .ped file is named X.ped, rename former files to X.ped.bed ...

13. to conduct the gwas, do

`./gemma -bfile X.ped -p PHENOTYPE.csv -lm -o OUTPUT`

.
.
.
.
Below is a catalog of all the errors in the process, not necessary for conducting the analysis

Steps (attempt 1):

1. fastqc applied to given fastq R1/R2 files. Primer cutting and quality filtering above score of 30 accomplished with

`cutadapt -q 30 -a PRIMER_OF_R1 -A PRIMER_OF_R2 -o OUTPUT1  -p OUTPUT2 INPUT1 INPUT2`

**process verified by the html output of a separate fastqc analysis**
i created a bash script which automated this process over 60 samples.

2. Downloaded reference genome of septoria musiva from [here](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000320565.1/)

    - in download options, only **☐ Genomic Sequence** was selected

3. Using [this](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format), a .dict file and .fai file were created from the reference file

`gatk CreateSequenceDictionary -R REFERENCE.fasta`

`samtools faidx REFERENCE.fasta`

4. Using [this](https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/BWA_tutorial.pdf) tutorial the reference genome was indexed with

`bwa index REFERENCE.FASTA`

5. Reads were combined into a SAM file using 

`bwa mem REFERENCE.fasta INPUT_R1.fastq INPUT_R2.fastq > OUTPUT.sam`

6. sam file converted to bam file using

`samtools view -h -b -S INPUT.sam > OUTPUT.bam`

7. filter bam file for only sequences that were mapped against reference genome

`samtools view -b -F 4 OUTPUT.bam > FILTERED.bam`

i created another bash script which automated steps 5-7 for the 60 samples given.

8. Using [this](https://gatk.broadinstitute.org/hc/en-us/articles/360035531892-GATK4-command-line-syntax) a vcf file creation was _attempted_ using

`gatk HaplotypeCaller -R REFERENCE.fasta -I FILTERED.bam -O OUTPUT.vcf`

but an error was returned:
![image](https://user-images.githubusercontent.com/108294550/178332364-c583023a-1213-458b-ba2d-736b13ad2f98.png)
the chief issue appears to be in the message:

"_java.lang.IllegalArgumentException: samples cannot be empty_"

.

Troubleshooting:

I googled the previous error and found [this](https://gatk.broadinstitute.org/hc/en-us/community/posts/360063062572-Not-getting-vcf-file) board. 

The answer suggested that it relates to an issue with the sam/bam file and encouraged [this](https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile) tutorial.

Being that the sam file creation was the first step, i ran ValidateSamFile on the sam file as

`gatk ValidateSamFile -I OUTPUT.sam`

and got a long list of errors and warnings:

![image](https://user-images.githubusercontent.com/108294550/178347769-f32532c0-88a3-4809-9db4-bc6dfc4062d0.png)
the key information is:

ERROR:MISSING_READ_GROUP	1

WARNING:QUALITY_NOT_STORED	5

WARNING:RECORD_MISSING_READ_GROUP	11589713

I tried this over many different samples (not just Nisk1 as shown) and got the same issue. Ricardo suspected that this is due to the oringinal fastq files not actually being septoria. Someone could have misnamed them, leading to all these issues. This would make sense of the line "WARNING:RECORD_MISSING_READ_GROUP	11589713" where it seems to suggest that the given genomes do not posses any of the read groups relative to the reference genome (which we are more certain of). 

To test this, i practiced assembling a genome into a fasta file from two read files (that were primer and quality filtered) after step 1. the following command accomplsihed this:
`spades.py -1 INPUT.R1.fastq -2 INPUT.R2.fastq -o OUTPUT.fasta`

.

we then applied bbsketch to the reference genome and the newly generated fasta. the bbsketch results confirmed that the fastq dataset given was indeed septoria. 

.

We were then back at square one. We tried googling the error and realized the error "ERROR:MISSING_READ_GROUP	1" simply means that its missing a read group, so we have to add them! 

I installed picard.jar to accomplish this and used the [following](https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) tutorial.

I guessed that adding the read group needed to be done after step 5 in the steps for attempt 1 because doing it after would force us to work in binary, which is not what the addorreplacegroups function works in. Also, ricardo hinted at this. so i took an output file from step 5 (sam file) and added a read group using:

`java -jar PICARD.jar AddOrReplaceReadGroups -I INPUT.sam -O OUTPUT.sam -RGLB lib1 -RGPL illumna -RGPU unit1 -RGSM 20`

I got a successful output file and continued with steps 6-8 to get the filtered bam file. 

I attempted step 8 and got the following error:

"A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
Please index all input files:"

And so i indexed the filtered bam file from the product of step 7 as:


`samtools index FILTERED.bam`


and reran step 8 and got a vcf file!

.

Troubleshooting (continued):

vcf errors!

The vcf file is not formatted correctly
. i suspect this is due to the java AddOrReplaceReadGroup function not actually adding a read group, but creating a separate read group file. This new separate file just contains the read group and not the rest of the genome, so it is lower in overall size. I then tried to merge the read group file (.sam format) with the original sam by using [merge](http://www.htslib.org/doc/samtools-merge.html)

`samtools merge -n ORIGINAL.sam -r READGROUP.sam -o CLEANED.sam`

before proceeding, i used ValidateSamFile on the Nisk1actual to see what potential sources of error could be

`gatk ValidateSamFile -I CLEANED.sam -MODE SUMMARY`

and got the following issues:
![image](https://user-images.githubusercontent.com/108294550/178591208-ba5f5032-3ac4-4b2c-b800-bd4236150262.png)


in short, nothing was fixed and a new error was added. the merge solution does not work. 

.

Ricardo suggested that i try a different reference genome file that he provided. I repeated steps 3-5, and used the addoreplacereadgroups function again because when he did it the file size wasn't reduced. After using the AddOrReplaceReadGroups function, the same issue was present (the unfiltered [original] is not the same size as the Nisk1.sam, which is the output)![image](https://user-images.githubusercontent.com/108294550/178601475-5c14a7d8-9367-40d9-b52f-dcee00ca3456.png)

The issue was that there was that the read files i was using to generate the initial sam file were quality and primer filtered. I attempted the same process starting from step 5 with the original, raw read files and got a successful sam file with a read group:![image](https://user-images.githubusercontent.com/108294550/178606345-b0709add-1942-4a7d-a7bb-f49de202beb9.png)

I proceeded to, at the request of ricardo, replace steps 6-7 with the following commands

`java -jar picard.jar SortSam -I ID.sam -O ID.bam -SORT_ORDER coordinate -CREATE_INDEX true`

`java -jar picard.jar MarkDuplicates -I ID.bam -O ID.mdup.bam -M Nisk1.mdup -ASSUME_SORT_ORDER coordinate -CREATE_INDEX true`

*-M creates a .mdup file to write out duplicates to.*

`java -jar picard.jar SortSam -I ID.mdup.bam -O ID.sorted.bam -SORT_ORDER coordinate -CREATE_INDEX true`

to my understanding, these commands accomplish the same functions, just with a different package. they also sorted the contents better. It also eliminates the need to index the .bam files.

after running the three commands, the final output file was ID.sorted.bam. this was used as the input for step 8.

STEPS UPDATED (attempt 2)

*this assumes you are given the fastq reads and the reference genome*

1. index the reference fasta

`gatk CreateSequenceDictionary -R REF.fasta`

`samtools faidx REF.fasta`

`bwa index REF.fasta`

*keep everything here in the Reference_indexing directory*

2. combine R1 & R2 of the reads into a SAM file

`bwa mem REF.fasta ID_R1.fastq ID_R2.fastq > ID.sam`

*produces .sam file*

3. add read group to the sam file

`java -jar DIRECTORY/picard.jar AddOrReplaceReadGroups -I ID.sam -O ID.rg.sam -RGID ID -RGLB ID -RGPL Illumna -RGPU ID -RGSM ID`

*produces rg.sam files*

4. sort the sam file, conver to bam file, and index output

`java -jar DIRECTORY/picard.jar SortSam -I ID.rg.sam -O ID.bam -SORT_ORDER coordinate -CREATE_INDEX true`

*produces .bam, .bai files*

5. Mark duplicates, write them out to .mdup, filter bam file, and index output

`java -jar DIRECTORY/picard.jar MarkDuplicates -I ID.bam -O ID.mdup.bam -M ID.mdup -ASSUME_SORT_ORDER coordinate -CREATE_INDEX true`

*produces .mdup.bam, .mdup.bai files*

6. Sort the marked up bam file and index output

`java -jar picard.jar SortSam -I ID.mdup.bam -O ID.sorted.bam -SORT_ORDER coordinate -CREATE_INDEX true`

*produces .sorted.bam, .sorted.bai files*

7. produce vcf files from the .sorted.bam files

`gatk HaplotypeCaller -R REF.fasta -I ID.sorted.bam -O ID.vcf -ERC GVCF -ploidy 1`

*produces .vcf files*

8. produce a combined vcf file using CombineGVCFs


`gatk CombineGVCFs -R REF.fasta -V ID.vcf -V ID2.vcf (...) -V IDN.vcf -O COMBINED.vcf`

a. if you have a lot of files, create a script that prints all of the file names delimited by " -V " 

*produced .vcf FILE*

9. produce a joint genotyping

`gatk GenotypeGVCFs -R REF.fasta -V COMBINED.vcf -O COMBINED.gvcf.vcf`

*produces .vcf FILE*

.

Summarizing the reference genomes

Septoria Musiva (Sphaerulina Musiva) https://www.ncbi.nlm.nih.gov/assembly/GCF_000320565.1/#/st

> Septoria
Genome size: 29,352,103 \
No of contigs: 458 \
N50: 167,873 \
No. of genes: 10,228 \ 

Populus trichocarpa https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Populus_trichocarpa/101/ https://www.ncbi.nlm.nih.gov/assembly/GCF_000002775.4

> Populus
Genome size: 434,289,848 \
No of contigs: 8,318 \
N50: 552,806 \
No. of genes: 37,272 \

.

mini-test:
**given the product of step 9 from attempt 2, i will try and conduct GWAS using GEMMA and will record the process here**

1. Generate .bed, .bim, .fam files

`plink --vcf X.vcf --allow-extra-chr`

2. Generate .map, .ped, .nosex, and .log files

`plink --vcf X.vcf --recode --out PREFIX --allow-extra-ch`

3. be sure to rename the .bed, .bim. and .fam files to contain the .ped file's name. i.e., if your .ped file is named X.ped, rename former files to X.ped.bed ...

4. to conduct the gwas, do

`./gemma -bfile X.ped -p PHENOTYPE.csv -lm -o OUTPUT`

