# -*- coding: utf-8 -*-
"""
================================
Helpers for phantom localization
================================
"""

import os.path as op
import numpy as np
import matplotlib.pyplot as plt

import mne

actual_pos = mne.dipole.get_phantom_dipoles('otaniemi')[0]

base_path = op.join(op.dirname(__file__), '..', '..', 'phantom')

maxfilter_options = [False, True, 'mne']
dipole_amplitudes = [20, 100, 200, 1000]
dipole_indices = [5, 6, 7, 8]


def plot_errors(errors, kind):
    # Visualize the result
    xs = np.arange(3)
    xticklabels = ('Raw', 'SSS', 'SSS$_{MNE}$')
    ylim = [0, 20]
    fig, axs = plt.subplots(5, 1, figsize=(4, 8))
    for ai, dipole_amplitude in enumerate([20, 100, 200, 1000]):
        ax = axs[ai]
        for di in range(len(dipole_indices)):
            ax.plot(xs, errors[:, ai, di], label='%d' % dipole_indices[di])
        ax.set(title='%d nAm' % dipole_amplitude, ylim=ylim, xticks=xs,
               ylabel='Error (mm)', xlim=[0, 2])
        ax.set(xticklabels=[''] * len(xs))
        if ai == 3:
            ax.set(xticklabels=xticklabels)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='upper center',
                      bbox_to_anchor=(0.5, -0.25), ncol=2)
        ax.grid(True)
    fig.tight_layout()
    axs[-1].set_visible(False)
    for ext in ('png', 'pdf'):
        plt.savefig(op.join('figures', ('phantom_errors_%s.' % kind) + ext))
    plt.show()


def get_fwd():
    # They all have approximately the same dev_head_t
    info = mne.io.read_info(base_path + '/1000nAm/dip05_1000nAm.fif')
    sphere = mne.make_sphere_model(r0=(0., 0., 0.), head_radius=0.080)
    src_fname = 'phantom-src.fif'
    if not op.isfile(src_fname):
        mne.setup_volume_source_space(
            subject=None, fname=None, pos=3.5, mri=None,
            sphere=(0.0, 0.0, 0.0, 80.0), bem=None, mindist=5.0,
            exclude=2.0).save(src_fname)
    src = mne.read_source_spaces(src_fname)
    fwd_fname = 'phantom-fwd.fif'
    if not op.isfile(fwd_fname):
        mne.write_forward_solution(fwd_fname, mne.make_forward_solution(
            info, trans=None, src=src, bem=sphere, eeg=False,
            meg=True))
    fwd = mne.read_forward_solution(fwd_fname)
    return src, fwd


def get_data(dipole_idx, dipole_amplitude, use_maxwell_filter, show=False):
    data_path = base_path + '/%dnAm/' % dipole_amplitude
    if use_maxwell_filter is True:
        fname = 'dip%02d_%dnAm_sss.fif' % (dipole_idx, dipole_amplitude)
    else:
        fname = 'dip%02d_%dnAm.fif' % (dipole_idx, dipole_amplitude)

    raw_fname = op.join(data_path, fname)
    raw = mne.io.read_raw_fif(raw_fname, preload=True, verbose='error')
    raw.info['bads'] = ['MEG2233', 'MEG2422', 'MEG0111']

    events = mne.find_events(raw, stim_channel='STI201')
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
        preload=True, add_eeg_ref=False)
    evoked = epochs.average()

    if show:
        evoked.plot(spatial_colors=True)
    if show:
        evoked.plot_joint()

    evoked.crop(0, None)
    sphere = mne.make_sphere_model(r0=(0., 0., 0.), head_radius=0.08)
    cov = mne.compute_covariance(epochs, tmax=-0.05)
    return epochs, evoked, cov, sphere
