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
    
    # print(results)
    print(type(results))
    print(len(results))
    print(len(results[0]))
    print(len(results[0][1]))
    print('---------------------------------------------------')
    # print(results[0][1])
    dfs = []
    dicts_along = []
    for result in results:
        dfs.append(result[0])
        dicts_along.append(result[1])
    
    print(len(dicts_along)) # 4, one list per job
    print(dicts_along[0])
    print(len(dicts_along[0])) # 10, one dict per event 
    print(dicts_along[0][0]) # one dict
    
    results_df = pd.concat(dfs, ignore_index=True)
    # print(results_df)
    
    # Save data
    s = stream.name 
    config_name = s[len(s) - s[::-1].find('/'):s.find('.yaml')]
    os.system('mkdir -p {}'.format(cfg['output_folder']))
    results_df.to_hdf('{}/{}_{}.hdf5'.format(cfg['output_folder'], cfg['file_name'], config_name), key='seed_{}'.format(cfg['rnd_state_seed']))
    
    for dict_ in dicts_along:
        for key in dict_.keys():
            df = pd.DataFrame({key: dict_[key]})
            df.to_hdf('{}/{}_{}.hdf5'.format(cfg['output_folder'], cfg['file_name'], config_name), key='seed_{}_key'.format(cfg['rnd_state_seed'], key))
    
    end = time.time()
    print('Duration: {} s'.format(end-start))


### ----------------- MAIN ----------------------
if __name__ == '__main__':
    
    main()
