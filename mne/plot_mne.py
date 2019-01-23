import pandas as pd

from phantom_helpers import plot_errors

methods = ['dipfit', 'lcmv', 'music']
# methods = ['dipfit', 'lcmv']

dfs = []
postfix = ""
# maxfilter = 'False'

# maxfilter = 'False'
maxfilter = 'True'
# maxfilter = 'mne'
# postfix = "_aston"
postfix = "_aston2"

for m in methods:
    data = pd.read_csv('results/phantom_errors_%s%s.csv' % (m, postfix))
    print(data)
    if 'Unnamed: 0' in data.columns:
        data = data.drop(['Unnamed: 0'], axis=1)
    data['maxfilter'] = data['maxfilter'].astype(str)
    data = data.query('maxfilter == "%s"' % maxfilter).copy()
    dfs.append(data)

errors = pd.concat(dfs, axis=0, ignore_index=True)

plot_errors(errors, "mne_%s" % maxfilter, postfix="",
            ylim=(0, 15), xkey='method', xlabel_mapping=None)
