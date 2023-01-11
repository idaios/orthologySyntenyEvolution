list=$1
l=$2
m=$3
for f in $(cat $list); do
    fafile=$(ls *.pep.all.fa | grep -i $f)
    echo $fafile
   ./seqfilter.py -i $fafile -l $l $m -U > ${f}_${l}$m.fa
done
