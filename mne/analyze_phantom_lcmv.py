# -*- coding: utf-8 -*-
"""
=============================
Elekta phantom data with LCMV
=============================
"""
# authors: Amit & Alex & Eric

from __future__ import print_function

import numpy as np
import mne

from phantom_helpers import (get_data, get_fwd, actual_pos, maxfilter_options,
                             dipole_amplitudes, dipole_indices, plot_errors)


errors = np.empty(
    (len(maxfilter_options), len(dipole_amplitudes), len(dipole_indices)))
src, fwd = get_fwd()

for ui, mf in enumerate(maxfilter_options):
    for ai, dipole_amplitude in enumerate(dipole_amplitudes):
        print(('Processing : %4d nAm : SSS=%s'
               % (dipole_amplitude, mf)).ljust(40), end='')
        for di, dipole_idx in enumerate(dipole_indices):
            epochs, evoked, cov, sphere = get_data(
                dipole_idx, dipole_amplitude, mf)
            # Do LCMV
            data_cov = mne.compute_covariance(epochs, tmin=0.)
            stc = mne.beamformer.lcmv(
                evoked, fwd, cov, data_cov, reg=0.01, pick_ori='max-power')
            idx_max = np.argmax(np.mean(stc.data, axis=1))
            vertno_max = stc.vertices[idx_max]
            pos = src[0]['rr'][vertno_max]
            errors[ui, ai, di] = 1e3 * np.linalg.norm(
                pos - actual_pos[dipole_idx - 1])
            if dipole_amplitude < 1000 and errors[ui, ai, di] > 20:
                raise RuntimeError
        print(np.round(errors[ui, ai], 1))

plot_errors(errors, 'lcmv')
