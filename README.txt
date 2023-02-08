I. General.
This software package performs hierarchical clustering of samples based on transcriptomic (or any other numerical) data and the subsequent extraction of signatures of the identified clusters. Clustering can be carried out with constraints (must-links and cannot-links). The package is based on the Constrained–Divisive and Tjala algorithms, written in C++ and Python, and run in the command line. Raw data can be stored in various file formats, such as the following:
- text formats of tabular data (.tsv, .csv, .dat, .counts),
- spreadsheet formats (.xlsx…),
- data exchange formats (.json…),
- statistical formats (.dta, .sav, .por),
- data storage formats (.h5, .parquet, .orc).
II. Requirements.
Visual Studio, MinGW, and Python have been installed.
III. Beginning of work.
In the first run, .cpp files must be compiled. One of the possible options is the following:
1. From the Start Menu, open the Visual Studio 2022 folder and run the x64 Native Tools Command Prompt Line.
2. Enter cl.
3. Go to the folder with the files to be compiled.
4. Compile files with the following command: cl /EHsc <file_name>.
IV. Clustering (clustering.py, make_clustering.cpp, compute_distances.cpp, get_hierarchy.cpp, build_dendrogram.py)
1. Go to the folder with the clustering code and request information about the ordinal numbers of distance metrics with the following command: python clustering.py. 
2. Parse the raw data with the following command: python clustering.py <distance metric number>. The command requests for paths to folders with raw data, generates a new folder with parsed data, and creates the empty folder clustering_results for future results. This folder is created within the same superfolder that contains the raw data.
1 – Euclidean distance
2 – Pearson distance
3 – Spearman distance
4 – Manhattan distance
N.B. The Pearson distance is the most suitable standard metric for the analysis of transcriptomic “big data.”
3. Perform clustering with the following command: make_clustering.exe clustering_results/ <distance metric number>. This command creates files with the results of each iteration in the clustering_results folder.
4. Create a hierarchical tree from the formed clusters using the following command: get_hierarchy.exe clustering_results. This command creates a file in the clustering_results folder. This file contains information on which clusters each sample belongs to.
5. Optionally, visualize the clustering results as a dendrogram with the following command: python build_dendrogram.py without parameters.
V. Signature extraction (count_variances.py, count_variances.cpp, tjala.py)
1. Move the above program files to the clustering_results folder
2. Change to the clustering_results folder
3. python count_variances.py total_data.csv
4. count_variances.exe
5. python tjala.py
