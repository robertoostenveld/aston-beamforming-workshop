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
from phantom_helpers import get_dataset

# base_path, postfix = get_dataset('aston')
base_path, postfix = get_dataset('')

maxfilter_options, dipole_amplitudes, dipole_indices, actual_pos, bads =\
    get_bench_params(base_path)

_, fwd = get_fwd(base_path)

columns = ['dipole_index', 'dipole_amplitude', 'maxfilter', 'error']


def run(da, di, mf):
    print(('Processing : %4d nAm (dip %d) : SSS=%s'
          % (da, di, mf)).ljust(42), end='')
    epochs, evoked, cov, sphere = get_data(
        base_path, di, da, mf, bads=bads)
    pos = mne.beamformer.rap_music(
        evoked, fwd, cov, n_dipoles=1, return_residual=False)[0].pos[0]
    error = 1e3 * np.linalg.norm(pos - actual_pos[di - 1])
    print(" Error=%s mm" % np.round(error, 1))
    return pd.DataFrame([(di, da, mf, error)], columns=columns)


parallel, prun, _ = parallel_func(run, n_jobs=4)
errors = parallel([prun(da, di, mf) for mf, da, di in
                   product(maxfilter_options, dipole_amplitudes,
                           dipole_indices)])
errors = pd.concat(errors, axis=0, ignore_index=True)

plot_errors(errors, 'music', postfix=postfix)
