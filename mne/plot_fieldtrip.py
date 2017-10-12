from scipy.io import loadmat
import pandas as pd

from phantom_helpers import plot_errors

# data = loadmat('../fieldtrip/report/20170913/table.mat')
# data = loadmat('../fieldtrip/report/tmp/table.mat')
data = loadmat('../fieldtrip/report/table_movement_tsss_mc.mat')

table = data['table']

try:
    methods = list(map(lambda x: str(x[0]), data['methods'][0]))
    dip = list(map(lambda x: str(x[0]), data['dip'][0]))
    amp = list(map(lambda x: str(x[0]), data['amp'][0]))
except Exception as e:
    methods = list(map(lambda x: str(x[0]), data['methods'][0]))
    dip = list(map(lambda x: str(x), data['dip'][0]))
    amp = list(map(lambda x: str(x), data['amp'][0]))

errors = []
columns = ['dipole_index', 'dipole_amplitude', 'maxfilter', 'error']
postfix = ""

for i, mi in enumerate(methods):
    method = methods[i]
    for j, di in enumerate(dip):
        if di.startswith("dip"):
            di = int(di[3:])
        else:
            di = int(di)
        for k, da in enumerate(amp):
            error = table[j, k, i]
            if da.endswith('nAm'):
                da = int(da[:-3])
            else:
                da = int(da)
            errors.append(pd.DataFrame([(di, da, method, error)], columns=columns))

errors = pd.concat(errors, axis=0, ignore_index=True)

plot_errors(errors, "fieldtrip1", postfix=postfix, ylim=(0, 15),
            xlabel_mapping=None)
