# Loop through 0 to 10 migration events (-m)
seq 0 1 10 | while read i
do
    echo "treemix -global -i input_treemix.frq.gz -m $i -root SA -k 500 -bootstrap -o out.$i  > treemix_$i.log
    
    Rscript /home/lws2/workspace/BDYT/zaqizaba/quntiyichuan/genek/Part5.treemix_psmc_smcpp/script/treemixVarianceExplained.R  out.$i > out.$i.explained
    
    Rscript /home/lws2/workspace/BDYT/zaqizaba/quntiyichuan/genek/Part5.treemix_psmc_smcpp/script/draw_treemix.R out.$i  pop.order  out.$i" > treemix_$i.sh
done