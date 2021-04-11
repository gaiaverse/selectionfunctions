import sys, pickle
sys.path.append('/home/andy/Documents/Research/software/')
import sqlutilpy, getdata, numpy as np, h5py

#catalogue = "gaia_edr3.gaia_source"
catalogue = "andy_everall.gaia3_rand100m"
query = f"""select floor(phot_g_mean_mag/0.1)*0.1 as magbin,
                    percentile_cont(0.5) within group (ORDER BY SQRT(phot_g_n_obs)*phot_g_mean_flux_error/phot_g_mean_flux) as med_gamp
                    from {catalogue}
                    where phot_g_mean_mag is not NULL
                    group by magbin"""
print(query)
med_amp = sqlutilpy.get(query, asDict=True, **getdata.sql_args)

print('Processing data')
order = np.argsort(med_amp['magbin'])
for key in med_amp.keys():
    med_amp[key] = med_amp[key][order]

savefile = 'median_gamp_edr3.h'
print(f'Saving data: {savefile}')
with h5py.File(savefile, 'w') as hf:
    for key in med_amp.keys():
        hf.create_dataset(key, data=med_amp[key])
