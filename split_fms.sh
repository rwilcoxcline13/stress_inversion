cat TestFocalMechanisms1.csv | awk '{ print $0 > "FM_"NR".txt"}'

for i in $(seq 41)
do
cat FM_1.txt FM_"$i".txt > FM_"$i"_final.csv
done

rm -rf FM_*.txt
