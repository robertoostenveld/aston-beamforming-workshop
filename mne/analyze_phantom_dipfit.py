# -*- coding: utf-8 -*-
"""
===================================
Elekta phantom data with dipole fit
===================================
"""
# authors: Amit & Alex & Eric

from __future__ import print_function

from itertools import product

import pandas as pd
import mne
from mne.parallel import parallel_func

from phantom_helpers import get_data, plot_errors, get_bench_params
from phantom_helpers import get_dataset, compute_error

base_path, postfix = get_dataset('aston')
# base_path, postfix = get_dataset('')

(maxfilter_options, dipole_amplitudes, dipole_indices, actual_pos,
    actual_ori, bads) = get_bench_params(base_path)


def run(da, di, mf):
    print(('Processing : %4d nAm (dip %d) : SSS=%s'
          % (da, di, mf)).ljust(42), end="")
    epochs, evoked, cov, sphere = get_data(
        base_path, di, da, mf, bads=bads)
    # Hack to only use gradiometers
    epochs.pick_types(meg='grad')
    evoked.pick_types(meg='grad')
    # Do the dipole fit
    t_peak = 66e-3  # only fit the largest peak
    evoked.crop(t_peak, t_peak)
    dip = mne.fit_dipole(evoked, cov, sphere)[0]
    pos = dip.pos[0]
    ori = dip.ori[0]
    gof = dip.gof[0]
    actual_params = dict(actual_pos=actual_pos[di - 1],
                         actual_ori=actual_ori[di - 1],
                         actual_amp=da / 2.)
    error = compute_error(di, pos, ori, 1e9 * dip.amplitude,
                          **actual_params)
    error['gof'] = gof
    error['maxfilter'] = mf
    return pd.DataFrame(error)

parallel, prun, _ = parallel_func(run, n_jobs=1)
errors = parallel([prun(da, di, mf) for mf, da, di in
                   product(maxfilter_options, dipole_amplitudes,
                           dipole_indices)])
errors = pd.concat(errors, axis=0, ignore_index=True)
errors['method'] = 'dipfit'

plot_errors(errors, 'dipfit', postfix=postfix)
