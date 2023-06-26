#! /bin/bash

cd /Users/pgutjahr/Documents/muon-deflection-paper/deflection/deflection_parametrizations

 

start=`date +%s`

for config in configs/*; do 
    echo ${config}
    python deflection_method_shower_dist.py ${config} &
done

wait
end=`date +%s`
runtime=$((end-start))
echo "This took $runtime sec"
