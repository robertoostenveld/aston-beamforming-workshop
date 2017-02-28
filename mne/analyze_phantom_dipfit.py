# -*- coding: utf-8 -*-
"""
==========================================
 Elekta phantom data with dipole fit
==========================================
"""
# authors: Amit & Alex & Eric

from __future__ import print_function

import numpy as np
import mne

from phantom_helpers import (get_data, plot_errors, mfs, dipole_amplitudes,
                             dipole_indices, actual_pos)

errors = np.empty((len(mfs), len(dipole_amplitudes), len(dipole_indices)))
for ui, mf in enumerate(mfs):
    for ai, dipole_amplitude in enumerate(dipole_amplitudes):
        print(('Processing : %4d nAm : SSS=%s'
               % (dipole_amplitude, mf)).ljust(40), end='')
        for di, dipole_idx in enumerate(dipole_indices):
            epochs, evoked, cov, sphere = get_data(
                dipole_idx, dipole_amplitude, mf)
            # Do the dipole fit
            t_peak = 66e-3  # only fit the largest peak
            evoked.crop(t_peak, t_peak)
            pos = mne.fit_dipole(evoked, cov, sphere)[0].pos[0]
            errors[ui, ai, di] = 1e3 * np.linalg.norm(
                pos - actual_pos[dipole_idx - 1])
        print(np.round(errors[ui, ai], 1))

plot_errors(errors, 'dipfit')
