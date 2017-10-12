# -*- coding: utf-8 -*-
"""
==============================
Elekta Phantom data with MUSIC
==============================
"""

# authors: Amit & Alex & Eric

from __future__ import print_function

from itertools import product

import numpy as np
import pandas as pd
import mne
from mne.parallel import parallel_func

from phantom_helpers import get_data, plot_errors, get_bench_params, get_fwd
from phantom_helpers import get_dataset, compute_error

base_path, postfix = get_dataset('aston')
# base_path, postfix = get_dataset('')

(maxfilter_options, dipole_amplitudes, dipole_indices, actual_pos,
    actual_ori, bads) = get_bench_params(base_path)

_, fwd = get_fwd(base_path)

columns = ['dipole_index', 'dipole_amplitude', 'maxfilter', 'error']


def run(da, di, mf):
    print(('Processing : %4d nAm (dip %d) : SSS=%s'
          % (da, di, mf)).ljust(42), end='')
    epochs, evoked, cov, sphere = get_data(
        base_path, di, da, mf, bads=bads)
    # Hack to only use gradiometers
    epochs.pick_types(meg='grad')
    evoked.pick_types(meg='grad')
    dip = mne.beamformer.rap_music(
        evoked, fwd, cov, n_dipoles=1, return_residual=False)[0]
    t_idx = np.argmax(dip.gof)
    pos = dip.pos[t_idx]
    ori = dip.ori[t_idx]
    gof = dip.gof[t_idx]
    amp = dip.amplitude[t_idx]
    actual_params = dict(actual_pos=actual_pos[di - 1],
                         actual_ori=actual_ori[di - 1],
                         actual_amp=da / 2.)
    error = compute_error(di, pos, ori, amp, **actual_params)
    error['gof'] = gof
    error['maxfilter'] = mf
    return pd.DataFrame(error, index=[0])


parallel, prun, _ = parallel_func(run, n_jobs=4)
errors = parallel([prun(da, di, mf) for mf, da, di in
                   product(maxfilter_options, dipole_amplitudes,
                           dipole_indices)])
errors = pd.concat(errors, axis=0, ignore_index=True)
errors['method'] = 'music'

plot_errors(errors, 'music', postfix=postfix)
