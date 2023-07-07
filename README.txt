I. General.
This software package performs hierarchical clustering of samples based on transcriptomic (or any other numerical) data and subsequent extraction of signatures from the identified clusters. Clustering can be performed with specified constraints (must-links and cannot-links). The package is based on the Constrained–Divisive and Tjala algorithms, written in C++ and Python, and operates in the command line. The input data can be stored in various file formats, such as:
- text formats for tabular data (.tsv, .csv, .dat, .counts)
- electronic spreadsheet formats (.xlsx...)
- data exchange formats (.json...)
- statistical formats (.dta, .sav, .por)
- data storage formats (.h5, .parquet, .orc).

II. Requirements.
Visual Studio, MinGW, and Python are installed.

III. Getting Started.
On the first run, the .cpp files need to be compiled. One possible option is as follows:
1. Open the Visual Studio 2022 folder in the Start Menu and launch the x64 Native Tools Command Prompt Line.
2. Enter cl
3. Navigate to the folder containing the files to be compiled.
4. Compile the files using a command like "cl /EHsc <file_name>."

IV. Clustering (clustering.py, make_clustering.cpp, compute_distances.cpp, get_hierarchy.cpp, build_dendrogram.py).
1. Navigate to the folder with the clustering code and request information about the order numbers of distance metrics using the command "python clustering.py."
2. Parse the input data using the command "python clustering.py 2" (where 2 is the order number of the distance metrics). The command prompts for paths to folders containing the input data, creates a new folder with parsed data, and creates an empty "clustering_results" folder for future results. This folder is created in the parent folder where the input data is located.
The order number of the distance metrics
1: Euclidean distance
2: Pearson metric
3: Spearman metric
4: Manhattan distance
N.B.: The Pearson metric is the most suitable standard metric for analyzing transcriptomic Big Data.
3. Perform clustering using a command like "make_clustering.exe clustering_results/ 2." This command creates result files for each iteration in the "clustering_results" folder.
4. Display the hierarchy of the created clusters using the command "get_hierarchy.exe clustering_results." This command creates a file in the "clustering_results" folder showing which clusters each sample belongs to.
5. Optionally, visualize the clustering results as a dendrogram using the command "python build_dendrogram.py" without parameters.

V. Signature Extraction (count_variances.py, count_variances.exe, extract_clusters.py, genes.txt, tjala.py, samples_extractor.py)
1. Move the aforementioned files to the "clustering_results" folder.
2. Navigate to the "clustering_results" folder.
3. Specify the list of genes whose transcription level is used to generate signatures using the commands "python count_variances.py None" or "python count_variances.py genes.txt." The "genes.txt" file specifies a subgroup of genes (e.g., surfaceome genes) used to create the signatures. The "genes.txt" file should be in UTF-8 format. The command "python count_variances.py None" operates the full list of genes.
4. Calculate the variances using the command "count_variances.exe total_data.csv."
5. Generate a list of subclusters and the samples comprising each cluster using the command "python extract_clusters.py. Here, the "exit" command terminates only the "extract_clusters.py" command, not the entire pipeline.
6. Generate signatures using the command "python tjala.py." You can specify cannot clusters, i.e., clusters that should not be covered by the signatures (e.g., non-cancer samples when generating signatures for cancer cells). The list of cannot clusters can be generated using the "extract_clusters.py" command. If a cannot cluster is specified, the command prompts for a threshold indicating the minimum proportion of samples in the cannot cluster that should not be covered by the signatures. The command generates a file named "final_rules.txt" containing the signatures. For each signature, it provides the proportion of samples covered by it in each cluster. For each term of the signature, it provides a step scale indicating the proportion of the cannot cluster that would be excluded with each change in its value. The scale is provided for values ranging from 10-fold changes to 100-fold changes with 10 equal steps. The scale is provided to assess the prospects of developing drugs based on the given signature.
7. Obtain a list of samples that comply with the selected signature using the command "python samples_extractor.py." It creates a file with the list of samples in the "clustering_results" folder.
