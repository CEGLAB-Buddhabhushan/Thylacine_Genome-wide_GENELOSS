for i in `ls -1`
do
wc -l $i/"$i".fa|awk -F' ' '$1==0 {print $2}'|sed 's/\.fa//g'; done|cut -f1 -d'/'> not_checked.lst
wc not_checked.lst
for i in `cat not_checked.lst`
do
rm -r $i;./../../Confirm_Loss.sh $i
done
