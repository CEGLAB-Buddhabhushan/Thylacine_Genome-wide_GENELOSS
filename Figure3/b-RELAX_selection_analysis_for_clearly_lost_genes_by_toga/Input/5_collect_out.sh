echo "transcript ID,Antechinus_flavipes,Dasyurus_viverrinus,Dromiciops_gliroides,Gracilinanus_agilis,Lagorchestes_hirsutus,Macropus_fuliginosus,Macropus_giganteus,Macrotis_lagotis,Monodelphis_domestica,Myrmecobius_fasciatus,Notamacropus_eugenii,Phalanger_gymnotis,Potorous_gilbertii,Pseudocheirus_occidentalis,Pseudochirops_corinnae,Pseudochirops_cupreus,Sminthopsis_crassicaudata,Thylacinus_cynocephalus,Trichosurus_vulpecula,Vombatus_ursinus"|sed 's/,/\t/g' > Gene_loss_history.tsv
for i in XM_*
do
Antechinus_flavipes=`grep "Antechinus_flavipes" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Dasyurus_viverrinus=`grep "Dasyurus_viverrinus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Dromiciops_gliroides=`grep "Dromiciops_gliroides" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Gracilinanus_agilis=`grep "Gracilinanus_agilis" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Lagorchestes_hirsutus=`grep "Lagorchestes_hirsutus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Macropus_fuliginosus=`grep "Macropus_fuliginosus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Macropus_giganteus=`grep "Macropus_giganteus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Macrotis_lagotis=`grep "Macrotis_lagotis" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Monodelphis_domestica=`grep "Monodelphis_domestica" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Myrmecobius_fasciatus=`grep "Myrmecobius_fasciatus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Notamacropus_eugenii=`grep "Notamacropus_eugenii" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Phalanger_gymnotis=`grep "Phalanger_gymnotis" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Potorous_gilbertii=`grep "Potorous_gilbertii" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Pseudocheirus_occidentalis=`grep "Pseudocheirus_occidentalis" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Pseudochirops_corinnae=`grep "Pseudochirops_corinnae" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Pseudochirops_cupreus=`grep "Pseudochirops_cupreus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Sminthopsis_crassicaudata=`grep "Sminthopsis_crassicaudata" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Thylacinus_cynocephalus=`grep "Thylacinus_cynocephalus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Trichosurus_vulpecula=`grep "Trichosurus_vulpecula" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
Vombatus_ursinus=`grep "Vombatus_ursinus" $i/sorted/"$i".orthology.loss_summ.tsv|cut -f3`
echo $i
echo -e "$i\t$Antechinus_flavipes\t$Dasyurus_viverrinus\t$Dromiciops_gliroides\t$Gracilinanus_agilis\t$Lagorchestes_hirsutus\t$Macropus_fuliginosus\t$Macropus_giganteus\t$Macrotis_lagotis\t$Monodelphis_domestica\t$Myrmecobius_fasciatus\t$Notamacropus_eugenii\t$Phalanger_gymnotis\t$Potorous_gilbertii\t$Pseudocheirus_occidentalis\t$Pseudochirops_corinnae\t$Pseudochirops_cupreus\t$Sminthopsis_crassicaudata\t$Thylacinus_cynocephalus\t$Trichosurus_vulpecula\t$Vombatus_ursinus" >> Gene_loss_history.tsv
done
