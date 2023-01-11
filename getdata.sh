list=$1

for f in $(cat $list);
do
    #echo $f
    rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_fasta/$f/pep/*all*.gz .
done

#gunzip *.gz
