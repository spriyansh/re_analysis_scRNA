# Load the required libraries
import scanpy as sc
import phate as ph
import os
import sys
import numpy as np
import re
import matplotlib.pyplot as plt

# Declare Paths
inPrefix = "setty_et_al/data/output/Azimuth/"
outPrefix = "setty_et_al/data/output/PHATE/"
os.makedirs(outPrefix, exist_ok=True)

# List all files in the directory
avail_files = os.listdir(inPrefix)

# Subset the AnnData Files
avail_files = [name for name in avail_files if re.search('.h5ad', name)]

# Divide
rep1_names = [name for name in avail_files if re.search('rep1', name)] # Male-35
rep2_names = [name for name in avail_files if re.search('rep2', name)] # Female-28
rep3_names = [name for name in avail_files if re.search('rep3', name)] # Female-19

# Set up Phate
phate_operator = ph.PHATE(n_jobs=4, n_components=3, verbose=False)

# Create a list for traversal
file_name_list = list([rep1_names, rep2_names, rep3_names])

# Iteratively Compute Phate
for i in file_name_list:
    
    # Typecast
    i = "".join(i)
    
    # Select the Name of the Replicate
    parts = i.split("_")
    j = parts[0] 
    
    if j == "rep1":
        individual = "1"
        age = "35"
        sex = "Male"
    elif j == "rep2":
        individual = "2"
        age = "28"
        sex = "Female"
    elif j == "rep3":
        individual = "3"
        age = "19"
        sex = "Female"
    
    # Read the input file with prefix
    file_name = inPrefix + i
    annDataI = sc.read_h5ad(file_name)
    
    # Compute Phate
    phateI = phate_operator.fit_transform(annDataI)
    
    # Add Phate dimensions to the Object
    annDataI.obsm['PHATE'] = phateI
    
    # Plotting
    sc.pl.embedding(annDataI, basis='PHATE', color='cell.type', show=False)

    # Add custom Title
    plt.title("Individual: " + individual + " Age: " + age + " Sex: "+ sex)
    
    # Save the Plot
    file_name = outPrefix + j + "_PHATE2.png"
    plt.savefig(file_name, dpi=1400, bbox_inches='tight')
    plt.close()
    
    # Write
    file_name = outPrefix + j + "_PHATE3.tsv"
    custom_header = "PHATE1\tPHATE2\tPHATE3"
    np.savetxt(file_name, phateI, delimiter="\t", fmt = '%e', header = custom_header)
    
    # Validate 
    print("Completed for " + j)
