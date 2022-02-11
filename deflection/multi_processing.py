'''
------> multi processing kann auch im jupyter notebook verwendet werden
--> kopiere multi_03_new und füge zusätzlich noch die funktion deflection per interaction ein um bei bedarf auch die ablenkung pro 
interaktion zu ermitteln 
'''

import numpy as np
import click 
import yaml
from multiprocessing import Pool 
import time
import os
import pandas as pd
from helper_functions_multi import muon_propagation_custom_multi



@click.command() 
@click.argument('cfg', type=click.Path(exists=True))

def main(cfg):

    # Read data
    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
    
    # Get random numbers
    rnd = np.random.RandomState(cfg['rnd_state_seed'])
    rnds = rnd.randint(0, 9999, size=cfg['n_jobs'])

    # Set configs for propagation
    args = []
    for i in range(cfg['n_jobs']):
        arg = {
            'E_i': cfg['E_i'], 
            'E_f': cfg['E_f'], 
            'deflection': cfg['deflection'],
            'n_events': cfg['n_events_per_job'],
            'e_cut': cfg['e_cut'],
            'v_cut': cfg['v_cut'], 
            'cont_rand': cfg['cont_rand'],
            'scattering_method': cfg['scattering_method'],
            'deflection_type': cfg['deflection_type'],
            'table_path': cfg['table_path'],
        } 
        arg['rnd_seed'] = rnds[i]
        args.append(arg) # [arg]


    print('Total jobs: ', cfg['n_jobs'])
    print('Total events', cfg['n_jobs'] * cfg['n_events_per_job'])
    
    # Start propagation
    start = time.time()
    pool = Pool(cfg['n_CPU']) # number CPU cores
    results = pool.map(muon_propagation_custom_multi, args)
    pool.close()
    pool.join()
    
    results_df = pd.concat(results, ignore_index=True)
    print(results_df)
    
    # Save data
    s = stream.name 
    config_name = s[len(s) - s[::-1].find('/'):s.find('.yaml')]
    os.system('mkdir -p {}'.format(cfg['output_folder']))
    results_df.to_hdf('{}/{}_{}.hdf5'.format(cfg['output_folder'], cfg['file_name'], config_name), key='seed_{}'.format(cfg['rnd_state_seed']))
    
    end = time.time()
    print('Duration: {} s'.format(end-start))


### ----------------- MAIN ----------------------
if __name__ == '__main__':
    
    main()
