for strength in 1000nAm 200nAm 100nAm 20nAm ; do
echo '# Strength' $strength

for location in dip05 dip06 dip07 dip08 ; do
echo
echo '## Location' $location

for metric in dics lcmv rv ; do
if [ ${metric} = dics ] ; then
echo '### DICS scan, log10(active/baseline)'
elif [ ${metric} = lcmv ] ; then
echo '### LCMV scan, log10(active/baseline)'
elif [ ${metric} = rv ] ; then
echo '### Residual variance scan, -log10()'
fi

# it looks like this ![](dip05_1000nAm_sss.fif_rv.png)
echo '![]('${location}_${strength}_sss.fif_${metric}.png')'

done
done
done
