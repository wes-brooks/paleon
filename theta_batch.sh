for ((i=1; i<=19 ; i++))
do
    Rscript code/findtheta-fulllocal.r $i &
    Rscript code/findtheta.r $i &
done
