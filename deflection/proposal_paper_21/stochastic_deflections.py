import numpy as np
from tqdm import tqdm
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os


### bremsstrahlung
def get_new_psi_brems(E, E_, rnd_state, is_degree=True, theta_star=1):
    epsilon = E - E_
    mu = 0.1057  # muon mass
    p = rnd_state.uniform(0, 1)
    r_max = np.minimum(1, E_/epsilon) * E * theta_star / mu
    a = p * r_max**2 / (1+r_max**2)
    r = np.sqrt(a/(1-a))
    theta_photon = mu / E * r
    theta_mu = epsilon / E_ * theta_photon
    
    if is_degree:
        return np.rad2deg(theta_mu)
    else:
        return theta_mu
    
def get_new_psi_brems_ginneken_Eq6(E, E_, Z, m=0.10566):
    nu = (E - E_) / (E - m)
    
    k_1 = 0.092 * E**(-1/3)
    k_2 = 0.052 / E * Z**(-1/4)
    k_3 = 0.22 * E**(-0.92)
    rms_theta = np.max([np.min([k_1 * np.sqrt(nu), k_2]), k_3 * nu])
    if nu <= 0.5:    
        return rms_theta
    if (nu > 0.5) & (rms_theta < 0.2):
        k_4 = 0.26 * E**(-0.91)
        n = 0.81 * E**(0.5) / (E**(0.5) + 1.8)
        rms_theta = k_4 * nu**(1+n) * (1-nu)**(-n)
        return rms_theta
    if (nu > 0.5) & (rms_theta >= 0.2):
        n = 0.81 * E**(0.5) / (E**(0.5) + 1.8)
        k_4 = 0.26 * E**(-0.91) 
        k_5 = k_4 * nu**(1+n) * (1-nu)**(0.5-n)
        rms_theta = k_5 * (1-nu)**(-0.5)
        return rms_theta
    else:
        print('no case choosen')



rnd_state = np.random.RandomState(33)

colors = cm.Set1(np.linspace(0, 1, 9))


nu = np.logspace(-6,np.log10(0.999), 15) # 500
v = 1 - nu

n_events = int(1e5)

m = 0.1057

E = 1e2 # 100 GeV
rms_thetas = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    angles = [get_new_psi_brems(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2)))
plt.plot(nu, rms_thetas, 'x', color=colors[0], alpha=1, label=r'$E = 100\,$GeV (Geant4)')


E = 1e3 # 1 TeV
rms_thetas = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    angles = [get_new_psi_brems(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2))) 
plt.plot(nu, rms_thetas, 'x', color=colors[2], alpha=1, label=r'$E = 1\,$TeV (Geant4)')


E = 1e5 # 100 TeV
rms_thetas = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    angles = [get_new_psi_brems(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2)))    
plt.plot(nu, rms_thetas, 'x', color=colors[1], alpha=1, label=r'$E = 100\,$TeV (Geant4)')


### more data points at high nu
nu = np.array([0.65, 0.85, 0.95, 0.99])
E = 1e2 # 100 GeV
rms_thetas = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    angles = [get_new_psi_brems(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2)))
plt.plot(nu, rms_thetas, 'x', color=colors[0], alpha=1)

E = 1e3 # 1 TeV
rms_thetas = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    angles = [get_new_psi_brems(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2))) 
plt.plot(nu, rms_thetas, 'x', color=colors[2], alpha=1)


E = 1e5 # 100 TeV
rms_thetas = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    angles = [get_new_psi_brems(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2)))    
plt.plot(nu, rms_thetas, 'x', color=colors[1], alpha=1)




nu = np.logspace(-6,np.log10(0.999), 500) 

E = 100 # GeV
m = 0.1057
E_ = E - nu * (E - m)
values = [get_new_psi_brems_ginneken_Eq6(E, E__, Z=1) for E__ in E_]
plt.plot(nu, np.rad2deg(values), '--', color=colors[0], alpha=0.6, label=r'$E = 100\,$GeV (v.G.)')

E = 1000 # GeV
E_ = E - nu * (E - m)
values = [get_new_psi_brems_ginneken_Eq6(E, E__, Z=1) for E__ in E_]
plt.plot(nu, np.rad2deg(values), '--', color=colors[2], alpha=0.6, label=r'$E = 1\,$TeV (v.G.)')

E = 1e5 # GeV
E_ = E - nu * (E - m)
values = [get_new_psi_brems_ginneken_Eq6(E, E__, Z=1) for E__ in E_]
plt.plot(nu, np.rad2deg(values), '--', color=colors[1], alpha=0.6, label=r'$E = 100\,$TeV (v.G.)')

plt.xlabel(r'$\nu$')
plt.ylabel(r'$\langle\theta^2\rangle^{1/2}\,/\,°$')
plt.xscale('log')
plt.yscale('log')
# plt.xlim(5e-3, 0.5)
# plt.ylim(5e-5, 1e2)
plt.legend(loc='lower right', ncol=1)
plt.tight_layout()
plt.savefig('paper_bremsstrahlung_rms_theta_100GeV_1TeV_100TeV_allNu_Z1.pdf', bbox_inches='tight')
plt.clf()


### ioniz
def get_new_psi_deltaE_NEW(E_mu, E_mu_prime, m_e, m_mu, is_degree=True):
    # m_e = 0.511
    # m_mu = 105.658
    assert E_mu > m_mu, 'incoming energy lower than muon mass'
    
    gamma = E_mu / m_mu
    epsilon_max = 2 * m_e * (gamma**2 - 1) / (1 + 2*gamma*m_e/m_mu + (m_e/m_mu)**2) 
    assert E_mu_prime > E_mu - epsilon_max
    
    p_mu = np.sqrt((E_mu + m_mu) * (E_mu - m_mu))
    p_mu_prime = np.sqrt((E_mu_prime + m_mu) * (E_mu_prime - m_mu))
    
    cos_theta = ((E_mu + m_e) * E_mu_prime - E_mu*m_e - m_mu**2) / (p_mu * p_mu_prime)
    theta_mu = np.arccos(cos_theta)
    if is_degree:
        return np.rad2deg(theta_mu)
    else:
        return theta_mu
    
def get_nu_max(E_mu):
    # Energy in MeV
    # v_max = 1 - nu_max = 1 - E_f / E_i
    m_e = 0.51099895 # const.m_e * const.c**2 / const.e / 1e6
    m_mu = 105.6583755 # 1.883531627e-28 * const.c**2 / const.e / 1e6

    gamma = E_mu / m_mu
    nu_max = 1 - 2 * m_e * (gamma**2 - 1) / (1 + 2 * gamma * m_e/m_mu + (m_e/m_mu)**2)  / E_mu
    return nu_max


nu = np.logspace(-5, 0, 10000000) # 10000000


m_mu = 105.6583755
m_e = 0.51099895
E = 1e5 # in MeV
E_ = E - nu * (E - m_mu)
gamma = E / m_mu
epsilon_max = 2 * m_e * (gamma**2 - 1) / (1 + 2*gamma*m_e/m_mu + (m_e/m_mu)**2) 
E_cut = E_[E_ >= (E - epsilon_max)]
# print(E - epsilon_max)
# print(epsilon_max / (E - m_mu))
thetas = [get_new_psi_deltaE_NEW(E, E__, m_e, m_mu) for E__ in tqdm(E_cut)]
plt.plot(nu[E_ >= (E - epsilon_max)], thetas, '-', color=colors[0], label=r'$E = 100\,$GeV')
# print('E_i = {} --> E_f > {}'.format(E, get_nu_max(E)*E))



E = 1e6 # in MeV
E_ = E - nu * (E - m_mu)
gamma = E / m_mu
epsilon_max = 2 * m_e * (gamma**2 - 1) / (1 + 2*gamma*m_e/m_mu + (m_e/m_mu)**2) 
E_cut = E_[E_ >= (E - epsilon_max)]
# print(E - epsilon_max)
# print(epsilon_max / (E - m_mu))
thetas = [get_new_psi_deltaE_NEW(E, E__, m_e, m_mu) for E__ in tqdm(E_cut)]
plt.plot(nu[E_ >= (E - epsilon_max)], thetas, '-', color=colors[2], label=r'$E = 1\,$TeV')
# print('E_i = {} --> E_f > {}'.format(E, get_nu_max(E)*E))



E = 1e8 # in MeV
E_ = E - nu * (E - m_mu)
gamma = E / m_mu
epsilon_max = 2 * m_e * (gamma**2 - 1) / (1 + 2*gamma*m_e/m_mu + (m_e/m_mu)**2) 
E_cut = E_[E_ >= (E - epsilon_max)]
# print(E - epsilon_max)
# print(epsilon_max / (E - m_mu))
thetas = [get_new_psi_deltaE_NEW(E, E__, m_e, m_mu) for E__ in tqdm(E_cut)]
plt.plot(nu[E_ >= (E - epsilon_max)], thetas, '-', color=colors[1], label=r'$E = 100\,$TeV')
# print('E_i = {} --> E_f > {}'.format(E, get_nu_max(E)*E))


plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\theta\,/\,$°')
plt.legend()
plt.tight_layout()
plt.savefig('paper_ioniz_thetas_100GeV_1TeV_100TeV.pdf')
plt.clf()


### electron pair production
def get_rms_theta_pairprod_exp(E, nu):
    norm = 1e3
    n = -1
    a = 8.9e-4
    b = 1.5e-5
    c = 0.032
    d = 1
    e = 0.1
    # m = 105.7 / norm
    m_e =   0.5110 / norm
    # E = E * const.e * 1e9 / const.c**2
    # nu = (E - E_) / (E - m)
    minimum = np.min([a * nu**(1/4) * (1 + b*E) + c * nu / (nu + d), e])
    theta = (2.3 + np.log(E)) * (1- nu)**n / E * (nu - 2 * m_e/E)**2 / nu**2 * minimum
    return theta


E = 1e2 # in GeV ausgedrückt
nu = np.logspace(-7, -0.001, 200)
nu = nu[nu>4 * 511e-6/(E)]
thetas = [get_rms_theta_pairprod_exp(E, nus) for nus in nu]
plt.plot(nu, np.rad2deg(thetas), '-', color=colors[0], label=r'$E = 100\,$GeV (v.G.)')

E = 1e3 # in GeV ausgedrückt
nu = np.logspace(-7, -0.001, 200)
nu = nu[nu>4 * 511e-6/(E)]
thetas = [get_rms_theta_pairprod_exp(E, nus) for nus in nu]
plt.plot(nu, np.rad2deg(thetas), '-', color=colors[2], label=r'$E = 1\,$TeV (v.G.)')


E = 1e5 # in GeV ausgedrückt
nu = np.logspace(-7, -0.001, 200)
nu = nu[nu>4 * 511e-6/(E)]
thetas = [get_rms_theta_pairprod_exp(E, nus) for nus in nu]
plt.plot(nu, np.rad2deg(thetas), '-', color=colors[1], label=r'$E = 100\,$TeV (v.G.)')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\langle\theta^2\rangle^{1/2}\,/\,°$')
plt.legend()
plt.tight_layout()
plt.savefig('paper_epair_rms_theta_100GeV_1TeV_100TeV.pdf')
plt.clf()


### photonuclear interaction 
def get_new_psi_nuclint(E, E_, rnd_state, is_degree=True, nu_min=None, nu_max=True):
    M = 0.9383 # Proton mass
    mu = 0.1057 # Muon mass
    if nu_min is not None:
        if E - E_ < nu_min:
            print('minimum energy transfer is 0.2 GeV')
            return 0
    if nu_max:
        nu_max = E - M / 2
        if E - E_ > nu_max:
            print('maximum energy transfer is (E - mass_nucleon/2)={}, but here: ', E - E_)
            return 0
    m_0=np.sqrt(0.4)
    p = rnd_state.uniform(0, 1)
    # nu = epsilon
    epsilon = E - E_
    y = epsilon / E
    t_max = 2 * M * epsilon
    t_min = (mu * y)**2 / (1 - y)
    t_1 = np.minimum(epsilon**2, m_0**2)
    t_p = (t_max * t_1) / ((t_max + t_1) * ((t_max * (t_min + t_1))\
                    / (t_min * (t_max + t_1)))**p - t_max)
    sin2 = (t_p - t_min) / (4 * (E * E_ - mu**2) - 2 * t_min)
    theta_mu = 2 * np.arcsin(np.sqrt(sin2))
    
    if is_degree:
        return np.rad2deg(theta_mu)
    else:
        return theta_mu
    
def get_rms_theta_nuclint_ginneken(E, E_, m=0.1057):
    nu = (E - E_) / (E - m)
    rms_theta = (0.39 / (E * (1 - nu))) * (np.sqrt(E) * nu * (1 - nu))**0.17 * (1 - 0.135/(E * nu))
    return rms_theta    

n_events = int(1e5) 

nu_min_Geant4 = 0.201

nu = np.logspace(-7, -5e-3, 25)
m = 0.1057

E = 1e2 # 100 GeV
E_prime = E - nu * (E - m)
rms_thetas = []
rms_thetas_err = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    if E - E_ < nu_min_Geant4:
        # print(e, E - E_)
        continue
    angles = [get_new_psi_nuclint(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2)))  
    rms_thetas_err.append(np.std(angles))
plt.plot(nu[(E - E_prime) >= nu_min_Geant4], rms_thetas, 'x', color=colors[0], label=r'$E = 100\,$GeV (Geant4)')
# plt.errorbar(nu[(E - E_prime) > nu_min_Geant4], rms_thetas, yerr=rms_thetas_err, color='gray', alpha=0.5, fmt='x', label=r'$E = 100\,$GeV (Geant4)')


E = 1e3 # 1 TeV
E_prime = E - nu * (E - m)
rms_thetas = []
rms_thetas_err = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    if E - E_ < nu_min_Geant4:
        # print(e, E - E_)
        continue
    angles = [get_new_psi_nuclint(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2)))
    rms_thetas_err.append(np.std(angles))
plt.plot(nu[(E - E_prime) >= nu_min_Geant4], rms_thetas, 'x', color=colors[2], label=r'$E = 1\,$TeV (Geant4)')
# plt.errorbar(nu[(E - E_prime) > nu_min_Geant4], rms_thetas, yerr=rms_thetas_err, color='slategrey', fmt='x', label=r'$E = 1\,$TeV (Geant4)')



E = 1e5 # 100 TeV
E_prime = E - nu * (E - m)
rms_thetas = []
rms_thetas_err = []
for e in tqdm(nu):
    rms = np.array([])
    E_ = E - e * (E - m)
    if E - E_ < nu_min_Geant4:
        # print(e, E - E_)
        continue
    angles = [get_new_psi_nuclint(E, E_, rnd_state, is_degree=True) for i in range(n_events)]
    rms = np.append(rms, angles)
    rms_thetas.append(np.sqrt(np.mean(rms**2))) 
    rms_thetas_err.append(np.std(angles))
plt.plot(nu[(E - E_prime) >= nu_min_Geant4], rms_thetas, 'x', color=colors[1], label=r'$E = 100\,$TeV (Geant4)')
# plt.errorbar(nu[(E - E_prime) > nu_min_Geant4], rms_thetas, yerr=rms_thetas_err, color='black', fmt='x', label=r'$E = 100\,$TeV (Geant4)')

nu = np.logspace(-7, -5e-3, 200)

E = 1e2 # GeV
E_ = E - nu * (E - m)
values = np.rad2deg(get_rms_theta_nuclint_ginneken(E, E_))
plt.plot(nu[values>0], values[values>0], '-', color=colors[0], alpha=1, label=r'$E = 100\,$GeV (v.G.)')

E = 1e3 # GeV
E_ = E - nu * (E - m)
values = np.rad2deg(get_rms_theta_nuclint_ginneken(E, E_))
plt.plot(nu[values>0], values[values>0], '-', color=colors[2], alpha=1, label=r'$E = 1\,$TeV (v.G.)')

E = 1e5 # GeV
E_ = E - nu * (E - m)
values = np.rad2deg(get_rms_theta_nuclint_ginneken(E, E_))
plt.plot(nu[values>0], values[values>0], '-', color=colors[1], alpha=1, label=r'$E = 100\,$TeV (v.G.)')


plt.yscale('log')
plt.xscale('log')
# plt.xlabel(r'$\epsilon$ (Geant4) and $\nu$ (v.G.)')
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\langle\theta^2\rangle^{1/2}\,/\,°$')
plt.ylim(1e-6, 100)
plt.legend(loc='upper left', ncol=2)
plt.tight_layout()
plt.savefig('paper_photonuclear_rms_theta_100GeV_1TeV_100TeV.pdf')