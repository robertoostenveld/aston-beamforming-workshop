# -*- coding: utf-8 -*-
"""
==========================================
 Elekta phantom data with dipole fit
==========================================

"""

# authors: Amit & Alex

# import matplotlib
# matplotlib.use('Agg')  # for testing don't use X server

import os.path as op
import warnings
import numpy as np
import matplotlib.pyplot as plt
import mne
from mne.io import read_raw_fif

base_path = op.join(op.dirname(__file__), '..', '..', 'phantom')
stim_channel = 'STI201'
badchan = ['MEG2233', 'MEG2422', 'MEG0111']


def evaluate_dipole_fit(dipole_amplitude, show=False,
                        dipole_indices=[5, 6, 7, 8], use_maxwell_filter=True):
    data_path = base_path + '/%dnAm/' % dipole_amplitude

    errors = []
    for dipole_idx in dipole_indices:
        if use_maxwell_filter is True:
            fname = 'dip%02d_%dnAm_sss.fif' % (dipole_idx, dipole_amplitude)
        else:
            fname = 'dip%02d_%dnAm.fif' % (dipole_idx, dipole_amplitude)

        print('Processing : %s' % fname)
        raw_fname = op.join(data_path, fname)
        with warnings.catch_warnings(record=True) as _:
            raw = read_raw_fif(raw_fname, add_eeg_ref=False, preload=True)
        raw.info['bads'] = badchan

        events = mne.find_events(raw, stim_channel=stim_channel)
        if show: raw.plot(events=events)

        if show: raw.plot_psd(tmax=np.inf, fmax=60, average=False)

        raw.fix_mag_coil_types()
        if use_maxwell_filter == 'mne':
            # Use Maxwell filtering from MNE
            raw = mne.preprocessing.maxwell_filter(raw, origin=(0., 0., 0.))
            if show: raw.plot(events=events)

        ###############################################################################
        # We know our phantom produces sinusoidal bursts below 25 Hz, so let's filter.
        raw.filter(None, 40., h_trans_bandwidth='auto', filter_length='auto',
                   phase='zero')
        if show: raw.plot(events=events)

        ###############################################################################
        # Now we epoch our data, average it
        tmin, tmax = -0.15, 0.1
        event_id = events[0, 2]
        epochs = mne.Epochs(raw, events, event_id, tmin, tmax, baseline=(None, -0.05),
                            preload=True, add_eeg_ref=False)
        evoked = epochs.average()

        if show: evoked.plot(spatial_colors=True)
        if show: evoked.plot_joint()

        ###############################################################################
        # Let's do some dipole fits. The phantom is properly modeled by a single-shell
        # sphere with origin (0., 0., 0.). We compute covariance, then do the fits.
        t_peak = 66e-3  # ~66 MS at largest peak
        sphere = mne.make_sphere_model(r0=(0., 0., 0.), head_radius=None)
        cov = mne.compute_covariance(epochs, tmax=-0.05)
        dip, _ = mne.fit_dipole(evoked.copy().crop(t_peak, t_peak), cov, sphere, n_jobs=1)

        this_pos = actual_pos[dipole_idx - 1]
        this_errors = 1e3 * np.linalg.norm(dip.pos[0] - this_pos)
        errors.append(this_errors)
        print("Error for dipole idx %d: %s" % (dipole_idx, errors[-1]))

    return errors


actual_pos = [[56.3, 0., 32.5], [47.6, 0., 27.5], [39., 0., 22.5], [30.3, 0., 17.5],
              [32.5, 0., 56.3], [27.5, 0., 47.6], [22.5, 0., 39.], [17.5, 0., 30.3],
              [0., -56.3, 32.5], [0., -47.6, 27.5], [0., -39., 22.5], [0., -30.3, 17.5],
              [0., -32.5, 56.3], [0., -27.5, 47.6], [0., -22.5, 39.], [0., -17.5, 30.3],
              [-56.3, 0., 32.5], [-47.6, 0., 27.5], [-39., 0., 22.5], [-30.3, 0., 17.5],
              [-32.5, 0., 56.3], [-27.5, 0., 47.6], [-22.5, 0., 39.], [-17.5, 0., 30.3],
              [0., 32.5, 56.3], [0., 27.5, 47.6], [0., 22.5, 39.], [0., 17.5, 30.3],
              [0., 56.3, 32.5], [0., 47.6, 27.5], [0., 39., 22.5], [0., 30.3, 17.5]]
actual_pos = np.array(actual_pos) / 1e3  # For otaniemi phantom (in meter)

for use_maxwell_filter in [True, False, 'mne']:
    plt.figure()

    dipole_indices = [5, 6, 7, 8]

    for dipole_amplitude in [20, 100, 200, 1000]:
        errors = evaluate_dipole_fit(dipole_amplitude, dipole_indices=dipole_indices,
                                     use_maxwell_filter=use_maxwell_filter)
        plt.plot(dipole_indices, errors, label='%dnAm' % dipole_amplitude)

    plt.title('Dipole fit errors%s' % (' (SSS)' if use_maxwell_filter else ''))
    plt.xlabel('Dipole index')
    plt.ylabel('Error (mm)')
    plt.ylim(0, 10.)
    plt.grid('on')
    plt.legend(loc='upper left')
    plt.show()
    suffix = '_mf' if use_maxwell_filter else ''
    suffix += '_mne' if use_maxwell_filter is 'mne' else ''
    plt.savefig('figures/phantom_errors_dipfit%s.pdf' % suffix)
    plt.savefig('figures/phantom_errors_dipfit%s.png' % suffix)
