for i in `seq 1 6` X; do grep "TRANSCRIPT" ../$i/TOGA_Tasmanian_devil_"$i"_Tasmanian_wolf_o2o/loss_summ_data.tsv|cut -f2,3|sed 's/_/\./2' >> TOGA.loss_summ_data.tsv; echo $i; done
