###step1###
## Chr01
#snpbinner \
#crosspoints -r 0.001 -l 43270923 \
#--input chr01.tsv --output chr01.csv
## Chr02
snpbinner \
crosspoints -r 0.001 -l 35937250 \
--input chr02.tsv --output chr02.csv
## Chr03
snpbinner \
crosspoints -r 0.001 -l 36413819 \
--input chr03.tsv --output chr03.csv
## Chr04
snpbinner \
crosspoints -r 0.001 -l 35502694 \
--input chr04.tsv --output chr04.csv
## Chr05
snpbinner \
crosspoints -r 0.001 -l 29958434 \
--input chr05.tsv --output chr05.csv
## Chr06
snpbinner \
crosspoints -r 0.001 -l 31248787 \
--input chr06.tsv --output chr06.csv
## Chr07
snpbinner \
crosspoints -r 0.001 -l 29697621 \
--input chr07.tsv --output chr07.csv
## Chr08
snpbinner \
crosspoints -r 0.001 -l 28443022 \
--input chr08.tsv --output chr08.csv
## Chr09
snpbinner \
crosspoints -r 0.001 -l 23012720 \
--input chr09.tsv --output chr09.csv
## Chr10
snpbinner \
crosspoints -r 0.001 -l 23207287 \
--input chr10.tsv --output chr10.csv
## Chr11
snpbinner \
crosspoints -r 0.001 -l 29021106 \
--input chr11.tsv --output chr11.csv
## Chr12
snpbinner \
crosspoints -r 0.001 -l 27531856 \
--input chr12.tsv --output chr12.csv

for i in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12
do
###step2###
snpbinner \
bins -l 5000 -n $i -i ${i}.csv -o ${i}.bin.csv 

###step3###
snpbinner \
visualize \ 
-s ${i}.tsv \
-c ${i}.csv \
-b ${i}.bin.csv \
-o visualize_${i}
done
