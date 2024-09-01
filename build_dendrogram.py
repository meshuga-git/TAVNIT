import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from os import listdir
from os.path import isdir, join
import matplotlib.pyplot as plt
import random
import warnings
warnings.filterwarnings('ignore')

# Reading folders from files
folders_upper = []
with open('cannot_links.txt', 'r') as f:
    folders_upper.extend(f.readlines())
with open('must_links.txt', 'r') as f:
    folders_upper.extend(f.readlines())
with open('no_links.txt', 'r') as f:
    folders_upper.extend(f.readlines())

folders_upper = [x.strip() for x in folders_upper]

folders_all = []
for fold in folders_upper:
    all_obj = listdir(fold)
    folders_cur = []
    for q in all_obj:
        if isdir(join(fold, q)):
            folders_cur.append(join(fold, q))
    if len(folders_cur) == 0:
        folders_cur.append(fold)
    folders_all.extend(folders_cur)

# Assign random colors to folders
number_of_colors = len(folders_all)
color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
         for i in range(number_of_colors)]

total_data = pd.read_csv('clustering_results/total_data.csv', header=None)
filenames = total_data[0]
total_data[0] = [0] + [np.nan] * (len(total_data) - 1)

# Mapping folders to colors
folders_colors = {folders_all[x]: color[x] for x in range(len(folders_all))}

clusters_colors = {}
for x in range(len(filenames)):
    wo_parsed = filenames[x].replace('_parsed_data', '')
    for k, v in folders_colors.items():
        if k in wo_parsed:
            clusters_colors[x] = v
            break

cur_cluster = len(total_data)
Z_list = []
cols_reversed = list(reversed(total_data.columns))
prev = total_data[cols_reversed[0]].values
cur_cluster = len(total_data)

d = {x: x for x in range(len(total_data))}
cols = {x: 1 for x in range(len(total_data))}
dist = 1

# Cluster merging loop
for x in cols_reversed[1:]:
    if dist % 1000 == 0:
        print(f'{dist} objects have been processed')

    cur = total_data[x].values
    fst, snd = -1, -1
    for k in list(zip(prev, cur)):
        if not np.isnan(k[1]):
            fst = int(k[1])
        if not np.isnan(k[0]) and np.isnan(k[1]):
            snd = int(k[0])
            break

    fst_real = d[fst]
    snd_real = d[snd]
    cols_new = cols[fst_real] + cols[snd_real]

    del d[snd]
    d[fst] = cur_cluster

    # Ensure proper color assignment
    clusters_colors[cur_cluster] = clusters_colors.get(
        fst_real, clusters_colors.get(snd_real, '#' + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)]))
    )

    cols[cur_cluster] = cols_new
    del cols[fst_real]
    del cols[snd_real]

    cur_cluster += 1
    Z_list.append([fst_real, snd_real, dist + 0.0, cols_new])
    dist += 1
    prev = cur

# Build and save dendrogram
print('building dendrogram...')
plt.figure(figsize=(150, 150))
xx = dendrogram(Z_list, link_color_func=lambda x: clusters_colors[x])
plt.savefig('dendrogram.png', format='png', bbox_inches='tight')
print("Done! Dendrogram saved to file \'dendrogram.png\'")
