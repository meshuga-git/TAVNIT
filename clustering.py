import os
import sys
import math
import pandas as pd

import pandas as pd
import openpyxl
import warnings
warnings.filterwarnings('ignore')

from os import listdir, mkdir
from os.path import isfile, join, isdir
from collections import defaultdict

def process_all_files(all_files):

    index = 0
    for k, v in all_files.items():
        for f in v:
            csv = None
            if '.counts' in f or '.tsv' in f:
                csv = pd.read_csv(f, header = None, sep = '\t')
            elif '.csv' in f:
                csv = pd.read_csv(f, header = None, sep = ',')
            elif '.xls' in f:
                csv = pd.read_excel(f)
            elif '.dat' in f:
                csv = pd.read_table(f, header = None, sep="\s+")
            elif '.json' in f:
                csv = pd.read_json(f, orient='records')
            elif '.dta' in f:
                csv = pd.read_stata(f)
            elif '.sav' in f:
                csv = pd.read_spss(f)
            elif '.por' in f:
                import pyreadstat
                csv = pyreadstat.read_por(f)
            elif '.h5' in f:
                csv = pd.read_hdf(f)
            elif '.parquet' in f:
                csv = pd.read_parquet(f)
            elif '.orc' in f:
                csv = pd.read_orc(f)

            if csv is None:
                print('Unsupported format: ' + f)
                continue

            calc_columns = list(csv.columns)[1:]
            if len(calc_columns) == 1:
                calc_columns = ['expression']
            csv.columns = ['gene_raw'] + calc_columns
            csv['is_gene'] = csv.gene_raw.apply(lambda x: 1 if True else 0)
            csv_filtered = csv[csv.is_gene == 1]
            csv_filtered['gene'] = csv_filtered.gene_raw.apply(lambda x: x if x.find('.') == -1 else x[:x.find('.')])
            for q in calc_columns:
                cur_table = csv_filtered[['gene', q]]
                index += 1
                if q != 'expression':
                    cur_table.columns = ['gene', 'expression']
                    cur_table.to_csv(join(folder_for_parsed_data, k + '_' + q + '.csv'), index=False)
                else:
                    cur_table.to_csv(join(folder_for_parsed_data, k + '_' + f[f.rfind('\\')+1:f.find('.')].replace(k, '').replace(' ', '') + '.csv'), index=False)
                print('Object number ' + str(index) + ' parsed to csv.')


def parse_data(mypath):
    folders = [f for f in listdir(mypath) if isdir(join(mypath, f))]
    files = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[0] != '.']
    
    all_files = defaultdict(list)

    for x in folders:
        files_list = listdir(join(mypath, x))
        for f in files_list:
            all_files[x].append(join(join(mypath, x), f))

    for x in files:
        all_files[x[:x.find('.')]].append(join(mypath, x))
        
    return all_files


if __name__ == '__main__':

    if len(sys.argv) < 2 or sys.argv[1] not in ['1', '2', '3', '4']: # distance metric not specified
        print('Distance metric required: 1 - euclidean, 2 - pearson, 3 - spearman, 4 - manhattan')
        sys.exit()

    folders_all = []

    cannot_links = ''

    print('Write paths to 2 folders with cannot-links (separated by space) or "None" (without quotes):\n')
    res = str(input())
    if res != 'None' and len(res.split(' ')) != 2:
        print('Incorrect input. Please try again.')
        sys.exit()

    if res != 'None':
        folders_all.append(res.split(' ')[0])
        folders_all.append(res.split(' ')[1])
        cannot_links = '\n'.join(folders_all)

    with open('cannot_links.txt', 'w') as f:
        f.write(cannot_links)

    must_links = ''

    print('Write paths to the folder with must-links or "None" (without quotes):\n')
    res = str(input())
    if len(res.split(' ')) != 1:
        print('Incorrect input. Please try again.')
        sys.exit()

    if res != 'None':
        folders_all.append(res)
        must_links = res

    with open('must_links.txt', 'w') as f:
        f.write(must_links)

    no_links = ''

    print('Write paths to the folder with no constraints or "None" (without quotes):\n')
    res = str(input())
    if len(res.split(' ')) != 1:
        print('Incorrect input. Please try again.')
        sys.exit()

    if res != 'None':
        folders_all.append(res)
        no_links = res

    with open('no_links.txt', 'w') as f:
        f.write(no_links)

    for folder in folders_all:
        # parsing all files from the specified folders with a nesting level no greater than 2
        all_files = parse_data(folder)
    
        # all individual objects will be loaded into this folder - csv files consisting of 2 columns: gene and expression
        folder_for_parsed_data = folder + '_parsed_data'
        if not isdir(folder_for_parsed_data):
            mkdir(folder_for_parsed_data)
    
        # parsing all data and writing the resulting object files into the folder_for_parsed_data
        process_all_files(all_files)

        # writing the path to all files from the folder_for_parsed_data into a text file
        parsed_files = sorted(listdir(folder_for_parsed_data))
        parsed_files_for_cpp = '\n'.join([join(folder_for_parsed_data, x) for x in parsed_files])
        with open(folder + '_parsed_files.txt', 'w') as f:
            f.write(parsed_files_for_cpp)

    # results of all clustering iterations will be loaded into this folder - text files containing (comma-separated) cluster numbers corresponding to each object
    folder_for_clustering_results = 'clustering_results'
    if not isdir(folder_for_clustering_results):
        mkdir(folder_for_clustering_results)

    # compiling a C++ program that writes data from all files into a matrix, calculates all distances, and performs clustering
    # os.system('"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\cl" make_clustering.cpp')

    # running this program, passing as parameters the list of file paths, and the distance metric
    # os.system("make_clustering.exe " + folder_for_clustering_results + '/ ' + sys.argv[1])
