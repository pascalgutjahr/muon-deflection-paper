#! /bin/bash

cd /Users/pgutjahr/Documents/muon-deflection-paper/deflection/deflection_parametrizations

 

start=`date +%s`

for config in $1/*; do 
    echo ${config}
    python deflection_method_shower_dist.py ${config} &
done

wait
end=`date +%s`
runtime=$((end-start))
echo "This took $runtime sec"

# time ./run_proposal.py "path_to_configs"