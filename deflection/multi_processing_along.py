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
from iteration_utilities import flatten
from helper_functions_multi import muon_propagation_custom_multi, muon_propagation_custom_multi_along



@click.command() 
@click.argument('cfg', type=click.Path(exists=True))

def main(cfg):

    # Read data
    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
        
        
    # Print settings    
    for key in cfg: 
        print("{}: {}".format(key, cfg[key]))    
    print('Total events', cfg['n_jobs'] * cfg['n_events_per_job'])   
        
    
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
        if 'get_data_along_track' in cfg and cfg['get_data_along_track'] == True:
            arg['get_data_along_track'] = True
        else: arg['get_data_along_track'] = False
        args.append(arg) # [arg]

    
    # Start propagation
    start = time.time()
    pool = Pool(cfg['n_CPU']) # number CPU cores
    results = pool.map(muon_propagation_custom_multi_along, args)
    pool.close()
    pool.join()
    

    

    dfs_final = []
    dicts_along = []
    for dict_ in results:
        dfs_final.append(dict_['final_data'])
        dicts_along.append(dict_['data_along'])
        
    # print(len(dicts_along)) # 4, number of jobs 
    # print(len(dicts_along[0])) # 10, number of events per job
    # print(len(dicts_along[0][0])) # 5, number of dict entries in dict_data_along_track
    # print(len(list(flatten(dicts_along)))) # 40, number of all events 
    dicts_along_flatten = list(flatten(dicts_along))
    
    dfs_final_concat = pd.concat(dfs_final, ignore_index=True)
    print(dfs_final_concat)
    
    # Save data
    print('Saving data...')
    s = stream.name 
    config_name = s[len(s) - s[::-1].find('/'):s.find('.yaml')]
    os.system('mkdir -p {}'.format(cfg['output_folder']))
    dfs_final_concat.to_hdf('{}/{}_{}.hdf5'.format(cfg['output_folder'], cfg['file_name'], config_name), key='seed_{}'.format(cfg['rnd_state_seed']))
    
    # Save data for data along track
    for dict_ in dicts_along_flatten:
        for key in dict_.keys():
            df = pd.DataFrame({key: dict_[key]})
            df.to_hdf('{}/{}_{}.hdf5'.format(cfg['output_folder'], cfg['file_name'], config_name), format='table', key='seed_{}_{}'.format(cfg['rnd_state_seed'], key), append=True, index=True)
    
    end = time.time()
    print('Duration: {} s'.format(end-start))


### ----------------- MAIN ----------------------
if __name__ == '__main__':
    
    main()
