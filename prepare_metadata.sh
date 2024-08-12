printf 'RF,RS,Id,NumId,Population\n' > temp/metadata.csv
for i in $(ls Salvelinus-alpinus-RAD-Seq/DataDelivery_2024-08-08_13-44-45_ngisthlm00980/files/P32003/ | grep -v "\." | grep "_" | head -n 1) ;
do for e in $(find Salvelinus-alpinus-RAD-Seq/DataDelivery_2024-08-08_13-44-45_ngisthlm00980/files/P32003/"$i" -type f -name "*.fastq.gz") ;
do ln -s "$e" temp/data/fastq
done
files=$(ls temp/data/fastq/*"$i"* | sed 's/^/\.\//g' | tr '\n' ',' | sed 's/,$/\n/g')
numid=$(printf "$i" | sed 's/$/\n/g' | sed 's/.*_//g')
printf "$files","$i","$numid",default,'\n' >> temp/metadata.csv
done


