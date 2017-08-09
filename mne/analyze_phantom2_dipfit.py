# -*- coding: utf-8 -*-
"""
==========================================
Aston Elekta phantom data with dipole fit
==========================================
"""
# authors: Amit & Alex & Eric

from __future__ import print_function

import os.path as op
import numpy as np
import matplotlib.pyplot as plt

import mne

# from phantom_helpers import (plot_errors, maxfilter_options,
#                              dipole_amplitudes, dipole_indices, actual_pos)

# from phantom_helpers import actual_pos

actual_pos = mne.dipole.get_phantom_dipoles('vectorview')[0]

base_path = op.join(op.dirname(__file__), '..', '..', 'phantom_aston')

maxfilter_options = [False, 'mne']
dipole_amplitudes = [20, 200, 1000]
dipole_indices = [5, 6, 7, 8, 9, 10, 11, 12]

# dipole_amplitudes = [20]
# dipole_indices = [5]


def plot_errors(errors, kind):
    # Visualize the result
    xticklabels = ('Raw', 'SSS$_{MNE}$')
    xs = np.arange(len(xticklabels))
    ylim = [0, 20]
    fig, axs = plt.subplots(len(dipole_amplitudes) + 1, 1, figsize=(4, 8))
    for ai, dipole_amplitude in enumerate(dipole_amplitudes):
        ax = axs[ai]
        for di in range(len(dipole_indices)):
            ax.plot(xs, errors[:, ai, di], label='%d' % dipole_indices[di])
        ax.set(title='%d nAm' % dipole_amplitude, ylim=ylim, xticks=xs,
               ylabel='Error (mm)', xlim=[0, len(xticklabels) - 1])
        ax.set(xticklabels=[''] * len(xs))
        if ai == len(dipole_amplitudes) - 1:
            ax.set(xticklabels=xticklabels)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='upper center',
                      bbox_to_anchor=(0.5, -0.25), ncol=2)
        ax.grid(True)
    fig.tight_layout()
    axs[-1].set_visible(False)
    for ext in ('png', 'pdf'):
        plt.savefig(op.join('figures', ('phantom2_errors_%s.' % kind) + ext))
    plt.show()


def get_fwd():
    # They all have approximately the same dev_head_t
    # info = mne.io.read_info(base_path + '/1000nAm/dip05_1000nAm.fif')
    info = mne.io.read_info(base_path +
                            '/Amp1000_IASoff/Amp1000_Dip5_IASoff.fif')
    sphere = mne.make_sphere_model(r0=(0., 0., 0.), head_radius=0.080)
    src_fname = 'phantom2-src.fif'
    if not op.isfile(src_fname):
        mne.setup_volume_source_space(
            subject=None, fname=None, pos=3.5, mri=None,
            sphere=(0.0, 0.0, 0.0, 80.0), bem=None, mindist=5.0,
            exclude=2.0).save(src_fname)
    src = mne.read_source_spaces(src_fname)
    fwd_fname = 'phantom2-fwd.fif'
    if not op.isfile(fwd_fname):
        mne.write_forward_solution(fwd_fname, mne.make_forward_solution(
            info, trans=None, src=src, bem=sphere, eeg=False,
            meg=True))
    fwd = mne.read_forward_solution(fwd_fname)
    return src, fwd


def get_data(dipole_idx, dipole_amplitude, use_maxwell_filter, show=False):
    data_path = base_path + '/Amp%d_IASoff/' % dipole_amplitude
    # if use_maxwell_filter is True:
    fname = 'Amp%d_Dip%d_IASoff.fif' % (dipole_amplitude, dipole_idx)
    # else:
    #     fname = 'dip%02d_%dnAm.fif' % (dipole_idx, dipole_amplitude)

    raw_fname = op.join(data_path, fname)
    raw = mne.io.read_raw_fif(raw_fname, preload=True, verbose='error')
    # raw.info['bads'] = ['MEG2233', 'MEG2422', 'MEG0111']
    raw.info['bads'] = ['MEG1323', 'MEG1133', 'MEG0613', 'MEG1032', 'MEG2313',
                        'MEG1133', 'MEG0613', 'MEG0111', 'MEG2423']

    raw.crop(20, None)
    events = mne.find_events(raw, stim_channel='SYS201')
    if show:
        raw.plot(events=events)
    if show:
        raw.plot_psd(tmax=np.inf, fmax=60, average=False)

    raw.fix_mag_coil_types()
    if use_maxwell_filter == 'mne':
        # Use Maxwell filtering from MNE
        raw = mne.preprocessing.maxwell_filter(raw, origin=(0., 0., 0.))
        if show:
            raw.plot(events=events)

    #######################################################################
    # We know our phantom produces sinusoidal bursts below 25 Hz, so let's
    # filter.
    raw.filter(None, 40., h_trans_bandwidth='auto', filter_length='auto',
               phase='zero')
    if show:
        raw.plot(events=events)

    #######################################################################
    # Now we epoch our data, average it
    tmin, tmax = -0.15, 0.1
    event_id = events[0, 2]
    epochs = mne.Epochs(
        raw, events, event_id, tmin, tmax, baseline=(None, -0.05),
        preload=True)
    evoked = epochs.average()

    if show:
        evoked.plot(spatial_colors=True)
    if show:
        evoked.plot_joint()

    evoked.crop(0, None)
    sphere = mne.make_sphere_model(r0=(0., 0., 0.), head_radius=0.08)
    cov = mne.compute_covariance(epochs, tmax=-0.05)
    print(fname)
    print(evoked.nave)
    return epochs, evoked, cov, sphere


errors = np.empty(
    (len(maxfilter_options), len(dipole_amplitudes), len(dipole_indices)))

for ui, mf in enumerate(maxfilter_options):
    for ai, dipole_amplitude in enumerate(dipole_amplitudes):
        print(('Processing : %4d nAm : SSS=%s'
               % (dipole_amplitude, mf)).ljust(40), end='')
        for di, dipole_idx in enumerate(dipole_indices):
            epochs, evoked, cov, sphere = get_data(
                dipole_idx, dipole_amplitude, mf, show=False)
            # Do the dipole fit
            t_peak = 66e-3  # only fit the largest peak
            evoked.crop(t_peak, t_peak)
            pos = mne.fit_dipole(evoked, cov, sphere)[0].pos[0]
            errors[ui, ai, di] = 1e3 * np.linalg.norm(
                pos - actual_pos[dipole_idx - 1])
        print(np.round(errors[ui, ai], 1))

plot_errors(errors, 'dipfit')
