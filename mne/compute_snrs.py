# -*- coding: utf-8 -*-
"""
===================================
Compute SNR for Elekta phantom data
===================================
"""
# authors: Amit & Alex & Eric

from __future__ import print_function

from itertools import product

import os.path as op
import pandas as pd
import mne
from mne.parallel import parallel_func

from phantom_helpers import get_data, get_bench_params, get_fwd
from phantom_helpers import get_dataset

base_path, postfix = get_dataset('aston')
# base_path, postfix = get_dataset('')

(maxfilter_options, dipole_amplitudes, dipole_indices, actual_pos,
    actual_ori, bads) = get_bench_params(base_path)

# dipole_indices = [11]
# dipole_amplitudes = [1000]

src, fwd = get_fwd(base_path)


def run(da, di, mf):
    print(('Processing : %4d nAm (dip %d) : SSS=%s'
          % (da, di, mf)).ljust(43), end='')
    epochs, evoked, cov, sphere = get_data(
        base_path, di, da, mf, bads=bads)

    evoked_pst = evoked.copy().crop(0.05, None)
    inverse_operator = mne.minimum_norm.make_inverse_operator(
        epochs.info, fwd, cov, loose=1, depth=0.199)
    snr, snr_est = mne.minimum_norm.estimate_snr(evoked_pst, inverse_operator)
    peak_ch, peak_time = evoked_pst.get_peak(ch_type='mag')
    tp = int((peak_time - evoked_pst.times[0]) * evoked_pst.info['sfreq'])
    snr_peak = snr[tp]

    print(" SNR: %2.3f" % snr_peak)
    snrs = dict(actual_dipole_amplitude=da,
                dipole_index=di)
    snrs['snr'] = snr_peak
    snrs['maxfilter'] = mf
    return pd.DataFrame(snrs, index=[0])

parallel, prun, _ = parallel_func(run, n_jobs=4)
snrs = parallel([prun(da, di, mf) for mf, da, di in
                 product(maxfilter_options, dipole_amplitudes,
                         dipole_indices)])
snrs = pd.concat(snrs, axis=0, ignore_index=True)

snrs.to_csv(op.join('results', 'phantom_snrs_%s.csv' % postfix), index=False)
