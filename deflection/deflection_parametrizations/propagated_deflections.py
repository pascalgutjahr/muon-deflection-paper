import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py

file_path = 'data/deflection_method_shower_dist/paper.hdf5'
f = h5py.File(file_path)
keys = [key for key in f.keys()]
f.close()

df_dict = {}
for key in keys:
    # print(key)
    df_dict[key] = pd.read_hdf(file_path, key=key)


bins = np.logspace(-3, 2, 30)

key = 'E1000_moliere_v_cut0_001_e_cutinf_MS_only_vG'
plt.hist(np.rad2deg(df_dict[key].deflection), bins=bins, histtype='step', label='MSM Only')
median = np.median(np.rad2deg(df_dict[key].deflection))
upper_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.95)
lower_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.05)
df_dict[key]['median'] = median 
df_dict[key]['upper_90'] = upper_90 
df_dict[key]['lower_90'] = lower_90 
print(key)
print('max deflection: ', np.max(np.rad2deg(df_dict[key].deflection)))
print(f'median deflection: {median}, upper 90%: {upper_90}, lower 90%: {lower_90}')

key = 'E1000_moliereinterpol_v_cut0_001_e_cutinf_MS_only_vG'
plt.hist(np.rad2deg(df_dict[key].deflection), bins=bins, histtype='step', label='MSMI Only')
median = np.median(np.rad2deg(df_dict[key].deflection))
upper_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.95)
lower_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.05)
df_dict[key]['median'] = median 
df_dict[key]['upper_90'] = upper_90 
df_dict[key]['lower_90'] = lower_90 
print(key)
print('max deflection: ', np.max(np.rad2deg(df_dict[key].deflection)))
print(f'median deflection: {median}, upper 90%: {upper_90}, lower 90%: {lower_90}')

key = 'E1000_none_v_cut0_001_e_cutinf_Stoch_only_vG'
plt.hist(np.rad2deg(df_dict[key].deflection), bins=bins, histtype='step', label='Stoch Only vG')
median = np.median(np.rad2deg(df_dict[key].deflection))
upper_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.95)
lower_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.05)
df_dict[key]['median'] = median 
df_dict[key]['upper_90'] = upper_90 
df_dict[key]['lower_90'] = lower_90 
print(key)
print('max deflection: ', np.max(np.rad2deg(df_dict[key].deflection)))
print(f'median deflection: {median}, upper 90%: {upper_90}, lower 90%: {lower_90}')

key = 'E1000_moliereinterpol_v_cut0_001_e_cutinf_vG'
plt.hist(np.rad2deg(df_dict[key].deflection), bins=bins, histtype='step', label='MSMI + Stoch vG')
median = np.median(np.rad2deg(df_dict[key].deflection))
upper_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.95)
lower_90 = np.quantile(np.rad2deg(df_dict[key].deflection), 0.05)
df_dict[key]['median'] = median 
df_dict[key]['upper_90'] = upper_90 
df_dict[key]['lower_90'] = lower_90 
print(key)
print('max deflection: ', np.max(np.rad2deg(df_dict[key].deflection)))
print(f'median deflection: {median}, upper 90%: {upper_90}, lower 90%: {lower_90}')


plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta\,/\,°$')
plt.legend() # bbox_to_anchor=(1,1))
plt.tight_layout()   
plt.savefig('paper_deflection_4_params.pdf')
plt.clf()



bins = np.linspace(0, 5000, 30)

key = 'E1000_moliere_v_cut0_001_e_cutinf_MS_only_vG'
plt.hist(np.clip(df_dict[key].lateral_disp, bins[0], bins[-1]), bins=bins, histtype='step', label='MSM Only')
ld_median = np.median(df_dict[key].lateral_disp)
ld_upper_90 = np.quantile(df_dict[key].lateral_disp, 0.95)
ld_lower_90 = np.quantile(df_dict[key].lateral_disp, 0.05)
df_dict[key]['ld_median'] = ld_median 
df_dict[key]['ld_upper_90'] = ld_upper_90 
df_dict[key]['ld_lower_90'] = ld_lower_90 
print(key)
print('max disp: ', np.max(df_dict[key].lateral_disp))
print(f'median disp: {ld_median}, upper 90%: {ld_upper_90}, lower 90%: {ld_lower_90}')

key = 'E1000_moliereinterpol_v_cut0_001_e_cutinf_MS_only_vG'
plt.hist(np.clip(df_dict[key].lateral_disp, bins[0], bins[-1]), bins=bins, histtype='step', label='MSMI Only')
ld_median = np.median(df_dict[key].lateral_disp)
ld_upper_90 = np.quantile(df_dict[key].lateral_disp, 0.95)
ld_lower_90 = np.quantile(df_dict[key].lateral_disp, 0.05)
df_dict[key]['ld_median'] = ld_median 
df_dict[key]['ld_upper_90'] = ld_upper_90 
df_dict[key]['ld_lower_90'] = ld_lower_90 
print(key)
print('max disp: ', np.max(df_dict[key].lateral_disp))
print(f'median disp: {ld_median}, upper 90%: {ld_upper_90}, lower 90%: {ld_lower_90}')

key = 'E1000_none_v_cut0_001_e_cutinf_Stoch_only_vG'
plt.hist(np.clip(df_dict[key].lateral_disp, bins[0], bins[-1]), bins=bins, histtype='step', label='Stoch Only vG')
ld_median = np.median(df_dict[key].lateral_disp)
ld_upper_90 = np.quantile(df_dict[key].lateral_disp, 0.95)
ld_lower_90 = np.quantile(df_dict[key].lateral_disp, 0.05)
df_dict[key]['ld_median'] = ld_median 
df_dict[key]['ld_upper_90'] = ld_upper_90 
df_dict[key]['ld_lower_90'] = ld_lower_90 
print(key)
print('max disp: ', np.max(df_dict[key].lateral_disp))
print(f'median disp: {ld_median}, upper 90%: {ld_upper_90}, lower 90%: {ld_lower_90}')

key = 'E1000_moliereinterpol_v_cut0_001_e_cutinf_vG'
plt.hist(np.clip(df_dict[key].lateral_disp, bins[0], bins[-1]), bins=bins, histtype='step', label='MSMI + Stoch vG')
ld_median = np.median(df_dict[key].lateral_disp)
ld_upper_90 = np.quantile(df_dict[key].lateral_disp, 0.95)
ld_lower_90 = np.quantile(df_dict[key].lateral_disp, 0.05)
df_dict[key]['ld_median'] = ld_median 
df_dict[key]['ld_upper_90'] = ld_upper_90 
df_dict[key]['ld_lower_90'] = ld_lower_90 
print(key)
print('max disp: ', np.max(df_dict[key].lateral_disp))
print(f'median disp: {ld_median}, upper 90%: {ld_upper_90}, lower 90%: {ld_lower_90}')

# plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'lateral displacement in cm')
plt.legend() # bbox_to_anchor=(1,1))
plt.tight_layout()
plt.savefig('paper_lateral_disp_4_params.pdf')