from scipy.io import loadmat
import pandas as pd

from phantom_helpers import plot_errors

data = loadmat('../fieldtrip/report/20170906/table.mat')

table = data['table']
methods = list(map(lambda x: str(x[0]), data['methods'][0]))
dip = list(map(lambda x: str(x[0]), data['dip'][0]))
amp = list(map(lambda x: str(x[0]), data['amp'][0]))

errors = []
columns = ['dipole_index', 'dipole_amplitude', 'maxfilter', 'error']
postfix = ""

for i, mi in enumerate(methods):
    method = methods[i]
    for j, di in enumerate(dip):
        di = int(di[3:])
        for k, da in enumerate(amp):
            error = table[j, k, i]
            da = int(da[:-3])
            errors.append(pd.DataFrame([(di, da, method, error)], columns=columns))

errors = pd.concat(errors, axis=0, ignore_index=True)

plot_errors(errors, "fieldtrip_" + method, postfix=postfix, ylim=(0, 25),
            xlabel_mapping=None)
