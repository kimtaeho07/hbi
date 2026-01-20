### Make output and back up directories 
ODIR=/oak/stanford/groups/tchelepi/kimtaeho/Fervo/hbi/hbi_gaussianSource/pSig1
if [ ! -e $ODIR ]; then
    mkdir $ODIR 
fi
ln -s $ODIR output
rm -r $ODIR/*

