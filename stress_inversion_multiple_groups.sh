cd ~/Desktop/stress_inversion_julia/Final_Stress_Inversion_072117

input1=(single_foc TestFocalMechanisms)

for i in "${input1[@]}"
do
  sed "s/foc_file/$i/" Inputfile_template.txt > Inputfile.txt
  julia stress_inversion_grouping_08292017.jl Inputfile.txt
done

cd ~/Desktop/FM_data

for i in "${input1[@]}"
do
  mkdir "$i"_data
  mv filename_"$i"* "$i"_data
done
