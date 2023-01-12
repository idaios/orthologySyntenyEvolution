echo -e "\t0\t5\t10\t50\t90\t95\t100"
for f in `ls original/*.fa`; do
    ./distrLengths.py -i $f;
done | sed 's/original\///g'| sed 's/\.[^\t]*//g'
