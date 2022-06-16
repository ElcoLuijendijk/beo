"""

"""

__author__ = 'Elco Luijendijk'

import os
import pandas as pd

# add the parameters that you wish to transfer to the reorganized output file here:
index = 'runtime_yr'
parameters_to_keep = ['exhumation_rate']

# select model result file
result_dir = 'model_output'
files = os.listdir(result_dir)
files = [os.path.join(result_dir, f) for f in files if f[-4:] == '.csv' and 'results' in f]
files.sort(key=os.path.getmtime)
files = files[::-1]

print('output files, from newest to oldest:')
for i, fn in enumerate(files):
    print(i, fn)

print('enter a number to select a file')

a = input()
fn = files[int(a)]

# read model data
df = pd.read_csv(fn)

# calculate half width
df['partial_reset_half_width'] = (df['x_max_partial_reset_surface'] - df['x_min_partial_reset_surface']) / 2.0
df['full_reset_half_width'] = (df['x_max_full_reset_surface'] - df['x_min_full_reset_surface']) / 2.0

# reshape table
parameters_to_keep = parameters_to_keep + ['partial_reset_half_width', 'full_reset_half_width']
dfp = df.pivot_table(index=index, columns='model_run', values=parameters_to_keep)

# save as new file
fn_new = fn.split('.')[0] + '_recombined.csv'
print('saving reorganized csv file as %s' % fn_new)
dfp.to_csv(fn_new)

print('done')




