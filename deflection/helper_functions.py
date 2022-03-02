import proposal as pp
from tqdm import tqdm
import numpy as np
import pandas as pd

### Muon propagation
def propagate_deflected_muons_custom(
    initial_energies, 
    minimum_energies, 
    inter_type=[
        pp.particle.Interaction_Type.ioniz, 
        pp.particle.Interaction_Type.brems, 
        pp.particle.Interaction_Type.photonuclear, 
        pp.particle.Interaction_Type.epair], 
    deflection='default', 
    deflection_type='m_scat+stochastic', 
    e_cut=500, 
    v_cut=0.05, 
    cont_rand=False, 
    scattering_method="highlandintegral", 
    beta_brems=1.0,
    beta_ioniz=1.0,
    beta_epair=1.0,
    beta_multiplescatter=1.0,
    beta_photonuclear=1.0,
    rnd_seed=1337, 
    initial_direction=[0, 0, 1], 
    table_path="/Users/pgutjahr/.cache/PROPOSAL"):

    '''Propagate muon tracks with deflection. Scaling of Bremsstrahlung opening angle can be done by beta.
    
    Parameters
    ----------
    initial_energies: list of energies
    minimum_energs: list of energies, lower propagation limit
    inter_type: list of interaction types for propagation/deflection, 
                default: ioniz, brems, nuclint, epairprod
    deflection_type: string, choose one:
            1. 'm_scat+stochastic'
            2. 'm_scat'
            3. 'stochastic'
    beta: scaling factor for chosen interaction type
    e_cut, v_cut, cont_rand: usual PROPOSAL energy cut settings
    initial_direction: list of initial direction (cartesian coordinates)
    table_path: string, path to interpolation tables
    '''
    pp.InterpolationSettings.tables_path = table_path   # version 7
    
    pp.RandomGenerator.get().set_seed(rnd_seed)
    args = {
            "particle_def": pp.particle.MuMinusDef(),
            "target": pp.medium.Ice(),
            "interpolate": True,
            "cuts": pp.EnergyCutSettings(e_cut, v_cut, cont_rand)
            }

    cross = pp.crosssection.make_std_crosssection(**args)
    multiple_scatter = pp.make_multiple_scattering(scattering_method, args["particle_def"], args["target"], cross, True)
    if deflection == 'default':
        print('Default deflection')
        stochastic_deflect = pp.make_default_stochastic_deflection(inter_type,
            args["particle_def"], args["target"])
    else:
        print('Costum deflection')
        stochastic_deflect = []
        for d in deflection:
            stochastic_deflect.append(pp.make_stochastic_deflection(d, 
            args["particle_def"], args["target"]))
        
    
    collection = pp.PropagationUtilityCollection()
    collection.displacement = pp.make_displacement(cross, True)
    collection.interaction = pp.make_interaction(cross, True)
    collection.time = pp.make_time(cross, args["particle_def"], True)
    collection.decay = pp.make_decay(cross, args["particle_def"], True)

    if deflection_type == 'stochastic':
        print('stochastic deflection')
        collection.scattering = pp.scattering.ScatteringMultiplier(
            stochastic_deflect, 
            [(pp.particle.Interaction_Type.brems, beta_brems), (pp.particle.Interaction_Type.ioniz, beta_ioniz), 
            (pp.particle.Interaction_Type.epair, beta_epair), (pp.particle.Interaction_Type.photonuclear, beta_photonuclear)])
    elif deflection_type == 'm_scat':
        print('multiple scattering')
        collection.scattering = pp.scattering.ScatteringMultiplier(multiple_scatter, beta_multiplescatter)
    elif deflection_type == 'm_scat+stochastic':
        print('multiple scattering and stochastic deflection')
        collection.scattering = pp.scattering.ScatteringMultiplier(
            multiple_scatter, 
            stochastic_deflect, # no list for default deflection!!!
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
    for E_i, E_min in zip(tqdm(initial_energies), minimum_energies):
        init_state.energy = E_i # initial energy in MeV
        track = prop.propagate(init_state, max_distance = 1e9, min_energy = E_min)
        tracks.append(track)
        
    return tracks


### get deflection angle
def get_angle_deviation(azimuth1, zenith1, azimuth2, zenith2, dtype='float128'):
    """Get opening angle of two vectors defined by (azimuth, zenith)
    Parameters
    ----------
    azimuth1 : 
        Azimuth of vector 1 in rad.
    zenith1 : 
        Zenith of vector 1 in rad.
    azimuth2 : 
        Azimuth of vector 2 in rad.
    zenith2 : 
        Zenith of vector 2 in rad.
    Returns
    -------
        The opening angle in rad between the vector 1 and 2.
        Same shape as input vectors.
    """
    cos_dist = (np.cos(azimuth1 - azimuth2, dtype=dtype) *
                np.sin(zenith1, dtype=dtype) * np.sin(zenith2, dtype=dtype) +
                np.cos(zenith1, dtype=dtype) * np.cos(zenith2, dtype=dtype))
    cos_dist = np.clip(cos_dist, -1., 1., dtype=dtype)
    return np.arccos(cos_dist, dtype=dtype)  


### name energy to save dataframe and plots
def energy_name(E):
    '''
        Parameters
        ----------
        E:
            energy in MeV
    '''
    e_unit = None
    if E / 1e12 >= 1:
        e_unit = '{}EeV'.format(int(E/1e12))
    elif E / 1e9 >= 1:
        e_unit = '{}PeV'.format(int(E/1e9))
    elif E / 1e6 >= 1:
        e_unit = '{}TeV'.format(int(E/1e6))
    elif E / 1e3 >= 1:
        e_unit = '{}GeV'.format(int(E/1e3))
    elif E >= 1:
        e_unit ='{}MeV'.format(int(E/1))
    elif E * 1e3 >= 1:
        e_unit = '{}KeV'.format(int(E/1e-3))
        
    assert e_unit != None, "energy too low"
    
    return e_unit


### get deflection per interaction
def get_zenith_deflections_along_track(tracks, param_name):
    '''
    Returns
    -------
        deflections in degree
    '''
    types = {
            'Interaction_Type.continuousenergyloss': 0,
            'Interaction_Type.epair': 1,
            'Interaction_Type.brems': 2,
            'Interaction_Type.photonuclear': 3,
            'Interaction_Type.ioniz': 4,
            'Interaction_Type.decay': 5,
    }
    
    defl_stoch = []
    defl_cont = []
    defl_type = []
    defl_angle_stoch = []
    defl_angle_cont = []
    
    for track in tqdm(tracks):
        zenith_last = track.track_directions()[0].spherical_coordinates[2]
        azimuth_last = track.track_directions()[0].spherical_coordinates[1]
        
        for typ, direction in zip(track.track_types()[1:], \
                                  track.track_directions()[1:]):
            zenith_new = direction.spherical_coordinates[2]
            azimuth_new = direction.spherical_coordinates[1]
            zenith_diff = abs(zenith_last - zenith_new)
            angle = get_angle_deviation(azimuth_last, zenith_last, azimuth_new, zenith_new)
            if str(typ) in ['Interaction_Type.epair',\
                            'Interaction_Type.brems',\
                            'Interaction_Type.photonuclear',\
                            'Interaction_Type.ioniz']:
                defl_stoch.append(zenith_diff)
                defl_type.append(types[str(typ)])
                defl_angle_stoch.append(angle)
            elif str(typ) == 'Interaction_Type.continuousenergyloss':
                defl_cont.append(zenith_diff)
                defl_type.append(types[str(typ)])
                defl_angle_cont.append(angle)
            # if particle decays, skip that interactionpoint
            # else:
                # print(typ)
                # defl_type.append(types[str(typ)])
            zenith_last = zenith_new
            azimuth_last = azimuth_new
    
    data_along_track = {
        '{}_along_defl_stoch'.format(param_name): np.rad2deg(defl_stoch),
        '{}_along_defl_cont'.format(param_name): np.rad2deg(defl_cont),
        '{}_along_defl_type'.format(param_name): np.array(defl_type),
        '{}_along_defl_angle_stoch'.format(param_name): np.rad2deg(defl_angle_stoch),
        '{}_along_defl_angle_cont'.format(param_name): np.rad2deg(defl_angle_cont),
        
    }
    return data_along_track


### save data per interaction
def save_data_along_dict(df_dir, hdf_name, dict_along):
    for key in dict_along:
        df = pd.DataFrame({key: dict_along[key]})
        df.to_hdf(df_dir + hdf_name, key=key)
        
        
### load data per interaction        
def load_data_along_dict(df_dir, hdf_name, param_name):
    dict_along = {}
    for key in ['defl_stoch', 'defl_cont', 'defl_type', 'defl_angle_stoch', 'defl_angle_cont']:
        df = pd.read_hdf(df_dir + hdf_name, key='{}_along_{}'.format(param_name, key))
        dict_along['{}_along_{}'.format(param_name, key)] = df['{}_along_{}'.format(param_name, key)].values
    return dict_along