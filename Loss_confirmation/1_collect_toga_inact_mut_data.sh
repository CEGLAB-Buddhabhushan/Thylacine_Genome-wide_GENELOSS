for i in `seq 1 6` X; do cat ../$i/TOGA_Tasmanian_devil_"$i"_Tasmanian_wolf_o2o/inact_mut_data.txt >> TOGA.inact_mut_data.txt; echo $i; done
sed -i 's/_/\./2' TOGA.inact_mut_data.txt
