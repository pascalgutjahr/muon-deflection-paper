import numpy as np
from tqdm import tqdm
from helper_functions import propagate_deflected_muons_custom, get_angle_deviation, energy_name, get_zenith_deflections_along_track, save_data_along_dict, load_data_along_dict

n_events = 1 # 5000 # 1
E_i = 1e9 # 1 PeV (corresponds to MeV)
E_f = 1e6 # 1 TeV (corresponds to MeV)
E_i_final = np.ones(n_events) * E_i
E_f_final = np.ones(n_events) * E_f

# table_path="/Users/pascalgutjahr/.cache/PROPOSAL" # lehrstuhl macbook
table_path="/net/nfshome/home/pgutjahr/.cache/PROPOSAL" # vollmond

hdf_name = 'df_{}_{}_{}events.hdf'.format(energy_name(E_i), energy_name(E_f), n_events)
# hdf_name

rnd_state = np.random.RandomState(1)
interpol_nodes = 200
scattering_method = "moliere"

### run data
param_name = 'default'
deflection = ['bremstsaiapproximation', 
              'ioniznaive',
              'photoborogpetrukhin',
              'epairginneken']

print("start propagation")
tracks_default = propagate_deflected_muons_custom(E_i_final, E_f_final, deflection=deflection, scattering_method=scattering_method, interpol_nodes=interpol_nodes, table_path=table_path)

print("start get angle")
deflection_default = [] 
for track in tqdm(tracks_default):
    d = get_angle_deviation(track.track_directions()[0].spherical_coordinates[1], 
                            track.track_directions()[0].spherical_coordinates[2], 
                            track.track_directions()[-1].spherical_coordinates[1], 
                            track.track_directions()[-1].spherical_coordinates[2])
    deflection_default.append(np.rad2deg(d))