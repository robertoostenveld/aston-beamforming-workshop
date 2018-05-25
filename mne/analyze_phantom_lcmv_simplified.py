# -*- coding: utf-8 -*-
"""
Created on Tue May 22 16:12:21 2018

@author: vlitvak
"""

import mne
import os
import numpy as np
import scipy.io

amp = 200
dip = 5
stim_channel = 'SYS201'
suffix = '_tsss'

raw_fname = '.\{}nAm\Amp{}_Dip{}_IASoff{}.fif'.format(amp, amp, dip, suffix)
#%%

os.chdir('C:\\home\\Data\\Aston\\newphantom')

raw = mne.io.read_raw_fif(raw_fname, preload=True, verbose='error')
raw.info['bads'] = ['MEG1133', 'MEG1323']

events = mne.find_events(raw, stim_channel=stim_channel)
raw.fix_mag_coil_types()

raw.filter(None, 40., h_trans_bandwidth='auto', filter_length='auto',
           phase='zero')

 

event_id = events[0, 2]

combined =  mne.Epochs(
raw, events, event_id, -0.1, 0.1, baseline=(-0.1, 0.1),
preload=True)


baseline = mne.Epochs(
raw, events, event_id, -0.1, 0, baseline=(-0.1, 0),
preload=True)


activation = mne.Epochs(
raw, events, event_id, 0, 0.1, baseline=(0, 0.1),
preload=True)

combined.pick_types(meg='grad')
activation.pick_types(meg='grad')
baseline.pick_types(meg='grad')
#%%
cov_combined = mne.compute_covariance(combined)
cov_act      = mne.compute_covariance(activation)
cov_bas      = mne.compute_covariance(baseline)


evoked_act = activation.average()
evoked_bas = baseline.average()

sphere = mne.make_sphere_model(r0=(0., 0., 0.), head_radius=0.1)

 
src = mne.setup_volume_source_space(
    subject=None, pos=10, mri=None,
    sphere=sphere, bem=None, mindist=0,
    exclude=0)

fwd = mne.make_forward_solution(
    raw.info, trans=None, src=src, bem=sphere, eeg=False,
    meg=True)
#%%

filters = mne.beamformer.make_lcmv(combined.info, fwd, cov_combined, reg=0.01, pick_ori='max-power', reduce_rank=True)

stc_act = mne.beamformer.apply_lcmv(evoked_act, filters)
stc_bas = mne.beamformer.apply_lcmv(evoked_bas, filters)

res = np.log10(np.mean(stc_act.data ** 2, axis=1))-np.log10(np.mean(stc_bas.data ** 2, axis=1))

res1 = np.mean(stc_act.data ** 2, axis=1);

scipy.io.savemat('frompy.mat', mdict={'res': res, 'pos': src[0]['rr'], 'in': src[0]['inuse'], 'res1': res1})

