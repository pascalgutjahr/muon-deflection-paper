'''
------> multi processing kann auch im jupyter notebook verwendet werden
--> kopiere multi_03_new und füge zusätzlich noch die funktion deflection per interaction ein um bei bedarf auch die ablenkung pro 
interaktion zu ermitteln 
'''

import numpy as np
from multiprocessing import Pool 
import time
import pandas as pd
from helper_functions_multi import muon_propagation_custom_multi





def start_propagation():
    
    E_i = 1e9 
    E_f = 1e6 
    
    deflection = ['tsaiapproximationbremsstrahlung', 
                  'naivionization',
                  'borogpetrukhinnuclearinteraction',
                  'kelnerpairproduction'] ### same as default
    e_cut = 500
    v_cut = 0.05
    cont_rand = False
    scattering_method="highlandintegral"
    deflection_type="m_scat+stochastic"
    
    # table_path="/Users/pgutjahr/.cache/PROPOSAL" # mac mini
    # table_path="/Users/pascalgutjahr/.cache/PROPOSAL" # lehrstuhl macbook
    table_path="/net/nfshome/home/pgutjahr/.cache/PROPOSAL" # vollmond

    file_name = "{}_{}".format(scattering_method, deflection_type) # sampled_energies
    output_folder = "data"
    
    n_CPU = 2
    n_jobs = 2
    n_events_per_job = 10
    total_events = n_jobs * n_events_per_job
    n_rnd = n_jobs
    rnd_state_seed = 3
    rnd = np.random.RandomState(rnd_state_seed)
    rnds = rnd.randint(0, 9999, size=n_rnd)

    args = []
    for i in range(n_rnd):
        arg = {
            'E_i': E_i, 
            'E_f': E_f, 
            'deflection': deflection,
            'n_events': n_events_per_job,
            'e_cut': e_cut,
            'v_cut': v_cut, 
            'cont_rand': cont_rand,
            'scattering_method': scattering_method,
            'deflection_type': deflection_type,
            'table_path': table_path,
        } 
        arg['rnd_seed'] = rnds[i]
        args.append(arg) # [arg]

    start = time.time()

    print('Total jobs: ', n_jobs)
    print('Total events', total_events)
    
    pool = Pool(n_CPU) # number CPU cores
    results = pool.map(muon_propagation_custom_multi, args)
    pool.close()
    pool.join()
    
    results_df = pd.concat(results, ignore_index=True)
    print(results_df)
    
    results_df.to_hdf('{}/{}.hdf5'.format(output_folder, file_name), key='seed_{}'.format(rnd_state_seed))
    
    end = time.time()
    print('Duration: {} s'.format(end-start))


### ----------------- MAIN ----------------------
if __name__ == '__main__':
    
    start_propagation()
