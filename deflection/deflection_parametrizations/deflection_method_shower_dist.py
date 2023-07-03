import proposal as pp
import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import sys 
sys.path.append('../')
from helper_functions_multi import get_angle_deviation
import os
import shutil
import click 
import yaml


@click.command()
@click.argument('cfg', type=click.Path(exists=True))

def main(cfg):
    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)


    data_dir = 'data/deflection_method_shower_dist/'
    os.system('mkdir -p {}'.format(data_dir))

    plot_dir = 'plots/deflection_method_shower_dist/'
    os.system('mkdir -p {}'.format(plot_dir))


    deflection = [
        'bremsginneken', # 'bremstsaiapproximation', 
        'ioniznaive',
        'photoginneken', # 'photoborogpetrukhin',
        'epairginneken'] ### same as default
        

    interpol_nodes = 100 # 200

    initial_direction = [0, 0, 1]

    beta_brems = 1.0
    beta_ioniz = 1.0
    beta_epair = 1.0
    beta_multiplescatter = 1.0
    beta_photonuclear = 1.0

       # remove old tables
    if cfg['remove_old_tables']:
        print(f"Removing {cfg['table_path']}")
        shutil.rmtree(cfg['table_path'])


    if 'table_path' in cfg:
        os.makedirs(cfg['table_path'], exist_ok=True)
        
        pp.InterpolationSettings.tables_path = cfg['table_path']      # version 7
        
    pp.InterpolationSettings.nodes_dndx_e = interpol_nodes
    pp.InterpolationSettings.nodes_dndx_v = interpol_nodes


    if cfg['e_cut'] == 'np.inf':
        e_cut = np.inf
    else:
        e_cut = cfg['e_cut']

    pp.RandomGenerator.get().set_seed(cfg['rnd_seed'])
    args = {
            "particle_def": pp.particle.MuMinusDef(),
            "target": pp.medium.Ice(),
            "interpolate": True,
            "cuts": pp.EnergyCutSettings(e_cut, cfg['v_cut'], cfg['cont_rand'])
            }

    cross = pp.crosssection.make_std_crosssection(**args)

    collection = pp.PropagationUtilityCollection()
    collection.displacement = pp.make_displacement(cross, True)
    collection.interaction = pp.make_interaction(cross, True)
    collection.time = pp.make_time(cross, args["particle_def"], True)
    collection.decay = pp.make_decay(cross, args["particle_def"], True)


    if 'MS_only' in cfg and cfg['MS_only']:
        multiple_scatter = pp.make_multiple_scattering(cfg['scattering_method'], args["particle_def"], args["target"], cross, True)
        collection.scattering = pp.scattering.ScatteringMultiplier(multiple_scatter, beta_multiplescatter)

    elif 'Stoch_only' in cfg and cfg['Stoch_only']:
        stochastic_deflect = []
        for d in deflection:
            stochastic_deflect.append(pp.make_stochastic_deflection(d, 
            args["particle_def"], args["target"]))

        collection.scattering = pp.scattering.ScatteringMultiplier(
            stochastic_deflect, 
            [(pp.particle.Interaction_Type.brems, beta_brems), (pp.particle.Interaction_Type.ioniz, beta_ioniz), 
            (pp.particle.Interaction_Type.epair, beta_epair), (pp.particle.Interaction_Type.photonuclear, beta_photonuclear)])

    else:
        multiple_scatter = pp.make_multiple_scattering(cfg['scattering_method'], args["particle_def"], args["target"], cross, True)
        stochastic_deflect = []
        for d in deflection:
            stochastic_deflect.append(pp.make_stochastic_deflection(d, 
            args["particle_def"], args["target"]))

        collection.scattering = pp.scattering.ScatteringMultiplier(
            multiple_scatter, 
            stochastic_deflect, 
            beta_multiplescatter, 
            [(pp.particle.Interaction_Type.brems, beta_brems), (pp.particle.Interaction_Type.ioniz, beta_ioniz), 
            (pp.particle.Interaction_Type.epair, beta_epair), (pp.particle.Interaction_Type.photonuclear, beta_photonuclear)])

    utility = pp.PropagationUtility(collection = collection)
    detector = pp.geometry.Sphere(pp.Cartesian3D(0,0,0), 1e20) # version 7
    density_distr = pp.density_distribution.density_homogeneous(args["target"].mass_density)


    prop = pp.Propagator(args["particle_def"], [(detector, utility, density_distr)])

    init_state = pp.particle.ParticleState()
    init_state.position = pp.Cartesian3D(0, 0, 0) # version 7
    init_state.direction = pp.Cartesian3D(initial_direction[0], initial_direction[1], initial_direction[2]) # version 7


    tracks = []
    for i in tqdm(range(cfg['n_events'])):
        init_state.energy = cfg['E_i'] # initial energy in MeV
        track = prop.propagate(init_state, max_distance = cfg['max_dist'], min_energy = cfg['E_f']) # max_dist=1e9
        tracks.append(track)


    tracks_zenith = np.array([])
    tracks_azimuth = np.array([])
    tracks_x = np.array([])
    tracks_y = np.array([])
    tracks_z = np.array([])
    tracks_pos_x = np.array([])
    tracks_pos_y = np.array([])
    tracks_pos_z = np.array([])
    tracks_propagated_distances = np.array([])
    tracks_deflection = np.array([])
    for track in tqdm(tracks):
        tracks_zenith = np.append(tracks_zenith, track.track_directions()[-1].spherical_coordinates[2])
        tracks_azimuth = np.append(tracks_azimuth, track.track_directions()[-1].spherical_coordinates[1])
        tracks_x = np.append(tracks_x, track.track_directions()[-1].x)
        tracks_y = np.append(tracks_y, track.track_directions()[-1].y)
        tracks_z = np.append(tracks_z, track.track_directions()[-1].z)
        tracks_pos_x = np.append(tracks_pos_x, track.track_positions()[-1].x)
        tracks_pos_y = np.append(tracks_pos_y, track.track_positions()[-1].y)
        tracks_pos_z = np.append(tracks_pos_z, track.track_positions()[-1].z)
        tracks_propagated_distances = np.append(tracks_propagated_distances, track.track_propagated_distances()[-1])
        tracks_deflection = np.append(tracks_deflection, get_angle_deviation(
            track.track_directions()[0].spherical_coordinates[1],
            track.track_directions()[0].spherical_coordinates[2],
            track.track_directions()[-1].spherical_coordinates[1],
            track.track_directions()[-1].spherical_coordinates[2],
            dtype='float64'
        ))


    def d_shower(x, y, in_meter=True):
        d = np.sqrt(x**2 + y**2)
        if in_meter:
            return d / 100 # distance in meter
        else:
            return  d 
    
    d = d_shower(tracks_pos_x, tracks_pos_y, in_meter=False)
        
    df = pd.DataFrame()
    df['x'] = tracks_x
    df['y'] = tracks_y
    df['z'] = tracks_z
    df['pos_x'] = tracks_pos_x
    df['pos_y'] = tracks_pos_y
    df['pos_z'] = tracks_pos_z
    df['zenith'] = tracks_zenith
    df['azimuth'] = tracks_azimuth
    df['propagated_distances'] = tracks_propagated_distances # in cm
    df['deflection'] = tracks_deflection # in rad
    df['lateral_disp'] = d # in cm

    key = 'E{}_{}_v_cut{}_e_cut{}'.format(int(cfg['E_i']/1e3), cfg['scattering_method'], cfg['v_cut'], e_cut).replace('.', '_') # E_i in GeV

    if 'MS_only' in cfg and cfg['MS_only']:
        key += '_MS_only'

    if 'Stoch_only' in cfg and cfg['Stoch_only']:
        key += '_Stoch_only'

    hdf_file = cfg['hdf_file']
    df.to_hdf(data_dir + f'{hdf_file}.hdf5', key=key)
    print(f'data save to {hdf_file} with the key {key}')

if __name__ == '__main__':
    main()
