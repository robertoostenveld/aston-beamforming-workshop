# -*- coding: utf-8 -*-
"""
==============================
Elekta phantom data with MUSIC
==============================
"""

# authors: Amit & Alex & Eric

from __future__ import print_function

import numpy as np
import mne

from phantom_helpers import (get_data, get_fwd, plot_errors, actual_pos,
                             maxfilter_options, dipole_indices,
                             dipole_amplitudes)


errors = np.empty(
    (len(maxfilter_options), len(dipole_amplitudes), len(dipole_indices)))

src, fwd = get_fwd()
for ui, use_maxwell_filter in enumerate(maxfilter_options):
    for ai, dipole_amplitude in enumerate(dipole_amplitudes):
        print(('Processing : %4d nAm : SSS=%s'
               % (dipole_amplitude, use_maxwell_filter)).ljust(40), end='')
        for di, dipole_idx in enumerate(dipole_indices):
            epochs, evoked, cov, sphere = \
                get_data(dipole_idx, dipole_amplitude, use_maxwell_filter)
            pos = mne.beamformer.rap_music(
                evoked, fwd, cov, n_dipoles=1, return_residual=False)[0].pos[0]
            errors[ui, ai, di] = 1e3 * np.linalg.norm(
                pos - actual_pos[dipole_idx - 1])
        print(np.round(errors[ui, ai], 1))

plot_errors(errors, 'music')
