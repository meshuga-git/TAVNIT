import pandas as pd
import numpy as np
import os

def zapoln(x):
    cur = None
    res = []
    for k in x:
        if not np.isnan(k):
            res.append(k)
            cur = k
        else:
            res.append(cur)
    return res

xx = pd.read_csv('total_data.csv', header=None)

for col in xx.columns[1:]:
    xx[col] = zapoln(xx[col])
    
d = dict()
d2 = set()

for k in xx.values:
    hier = set()
    for q in k[1:]:
        hier.add(int(q))
    d[k[0].split('/')[-1]] = ' '.join([str(len(hier))] + [str(x) for x in sorted(list(hier))])
    d2.add(k[0])
    
all_genes = set()

for k in d2:
    aa = pd.read_csv(k)
    all_genes.update(aa.gene)
    
merged = pd.DataFrame(sorted(list(all_genes)), columns = ['gene'])

for k in d2:
    aa = pd.read_csv(k)
    aa.columns = ['gene', k.split('/')[-1]]
    merged = merged.merge(aa, how='left')
    
merged.fillna(0, inplace=True)

genes_needed = []

for x in merged.values:
    if sum(x[1:]) != 0:
        genes_needed.append(x[0])
        
merged_needed = pd.DataFrame(genes_needed, columns=['gene']).merge(merged)

with open('hierarchy_new.txt', 'w') as f:
    f.write('\n'.join(pd.DataFrame(merged_needed.columns[1:]).merge(pd.DataFrame(d.items()))[1].values))

with open('columns_new.txt', 'w') as f:
    f.write('\n'.join(merged_needed.columns[1:]))
    
with open('genes_new.txt', 'w') as f:
    f.write('\n'.join(merged_needed['gene']))
    
dropped = merged_needed.drop('gene', axis=1)
dropped = dropped.astype(int)
dropped.to_csv('values_new.txt', header=None, index=False, sep=' ')

os.system("g++ count_variances.cpp")
os.system("./a.out")