import sys, pickle
sys.path.append('/home/andy/Documents/Research/software/')
import sqlutilpy, getdata, numpy as np, h5py

catalogue = "gaia_edr3.gaia_source"#"andy_everall.gaia3_rand100m"
query = f"""select floor(phot_g_mean_mag/0.1)*0.1 as magbin,
                    sum(phot_g_n_obs) as k,
                    sum(matched_transits) as n,
                    sum(astrometric_matched_transits) as n_ast,
                    count(*)
            from {catalogue}
            where phot_g_mean_mag is not NULL
            group by magbin"""
print(query)
obs_counts = sqlutilpy.get(query, asDict=True, **getdata.sql_args)

order = np.argsort(obs_counts['magbin'])
for key in obs_counts.keys():
    obs_counts[key] = obs_counts[key][order]

obs_counts['mean_eff'] = (obs_counts['k']+1)/(obs_counts['n']*9.+2)

savefile = 'expected_gobs_efficiency_edr3.h'
print(f'Saving data: {savefile}')
with h5py.File(savefile, 'w') as hf:
    for key in obs_counts.keys():
        hf.create_dataset(key, data=obs_counts[key])
