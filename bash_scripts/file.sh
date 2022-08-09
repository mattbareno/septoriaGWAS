ls ../data/Septoria/reads/*.gz > ../output/intermediate/allfiles.txt
#cuts the files by unique identifiers
cut -d "/" -f 5 ../output/intermediate/allfiles.txt > ../output/intermediate/int.txt
rm -r ../output/intermediate/allfiles.txt
cut -d "_" -f 1-2 ../output/intermediate/int.txt > ../output/intermediate/int2.txt
rm -r ../output/intermediate/int.txt
cut -d "-" -f 7-9 ../output/intermediate/int2.txt | uniq -d  > ../output/intermediate/identifier.txt
rm -r ../output/intermediate/int2.txt

