import numpy as np
import pandas as pd
import proposal as pp


### ------------- DEFINE FUNCTIONS ------------



### Muon propagation
def propagate_deflected_muons_custom_settings_multi( 
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
        # print('Default deflection')
        stochastic_deflect = pp.make_default_stochastic_deflection(inter_type,
            args["particle_def"], args["target"])
    else:
        # print('Costum deflection')
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
        # print('stochastic deflection')
        collection.scattering = pp.scattering.ScatteringMultiplier(
            stochastic_deflect, 
            [(pp.particle.Interaction_Type.brems, beta_brems), (pp.particle.Interaction_Type.ioniz, beta_ioniz), 
            (pp.particle.Interaction_Type.epair, beta_epair), (pp.particle.Interaction_Type.photonuclear, beta_photonuclear)])
    elif deflection_type == 'm_scat':
        # print('multiple scattering')
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
    
    
    return init_state, prop

def muon_propagation_custom_multi(args):

    # args = args[0]
    
    init_state, prop = propagate_deflected_muons_custom_settings_multi(deflection=args['deflection'], rnd_seed=args['rnd_seed'], table_path=args['table_path'], e_cut=args['e_cut'], v_cut=args['v_cut'], cont_rand=args['cont_rand'], scattering_method=args['scattering_method'], deflection_type=args['deflection_type'])
    
    E_i = args['E_i']
    E_min = args['E_f']
    n_events = args['n_events']
    
    E_i_l = []
    E_f_track_l = []
    distance_l = []
    deflection_l = []
    x_f_l = []
    y_f_l = []
    z_f_l = []
    for i in range(n_events):
        init_state.energy = E_i # initial energy in MeV
        track = prop.propagate(init_state, max_distance = 1e9, min_energy = E_min)
        # Prepare data
        E_f_track = track.track_energies()[-1] 
        distance = track.track_propagated_distances()[-1]
        deflection = get_angle_deviation(track.track_directions()[0].spherical_coordinates[1], track.track_directions()[0].spherical_coordinates[2], track.track_directions()[-1].spherical_coordinates[1], track.track_directions()[-1].spherical_coordinates[2])
        E_i_l.append(E_i)
        E_f_track_l.append(E_f_track)
        distance_l.append(distance)
        deflection_l.append(deflection)
        x_f_l.append(track.track_positions()[-1].x)
        y_f_l.append(track.track_positions()[-1].y)
        z_f_l.append(track.track_positions()[-1].z)
    
    # Save data
    df = pd.DataFrame()
    df['E_i'] = np.array(E_i_l) / 1e3 # in GeV
    df['E_f'] = np.array(E_f_track_l) / 1e3 # in GeV
    df['distances'] = np.array(distance_l) / 100 # in meter
    df['deflection'] = np.rad2deg(deflection_l) # in degree
    df['x_dir_i'] = init_state.direction.x * np.ones(n_events)
    df['y_dir_i'] = init_state.direction.y * np.ones(n_events)
    df['z_dir_i'] = init_state.direction.z * np.ones(n_events)
    df['x_i'] = init_state.position.x * np.ones(n_events) # position in cm
    df['y_i'] = init_state.position.y * np.ones(n_events)
    df['z_i'] = init_state.position.z * np.ones(n_events)
    df['x_f'] = x_f_l
    df['y_f'] = y_f_l
    df['z_f'] = z_f_l
        
    print('-- job done --')

    return df





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


