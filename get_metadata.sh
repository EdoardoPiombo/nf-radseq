for i in $(ls data/fastq/ | sed 's/_L006_R.*//g' | sort -u) ; do
readlink -f data/fastq/"$i"_L006_R1_001.fastq.gz >> column1.txt
readlink -f data/fastq/"$i"_L006_R2_001.fastq.gz >> column2.txt
grep "$i" doc/metadata.txt | cut -f 2 >> column3.txt
grep "$i" doc/metadata.txt | cut -f 3 >> column4.txt
grep "$i" doc/metadata.txt | cut -f 4 >> column5.txt
done

paste column1.txt column2.txt column3.txt column4.txt column5.txt | sed 's/\t/,/g' > data/metadata.csv

