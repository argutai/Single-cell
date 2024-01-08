import scanpy as sc
import numpy as np
import sys

path = sys.argv[1]  #N22_WTA-1.h5ad

### --- Read --- ###
OBJ = sc.read_h5ad(path)

OBJ.var_names_make_unique()

# Filter by reads and UMIs
sc.pp.filter_cells(OBJ, min_genes=50) # 100
sc.pp.filter_genes(OBJ, min_cells=3) # 3

# Annotate mitochondrial genes, and add a bunch of columns
OBJ.var['mt'] = OBJ.var_names.str.startswith('mt-') # print(sum(OBJ.var['mt']))
sc.pp.calculate_qc_metrics(OBJ, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Crop by upper and lower lims of n_genes_by_counts
upper_lim = np.quantile(OBJ.obs.n_genes_by_counts.values, 0.98) # 98
lower_lim = np.quantile(OBJ.obs.n_genes_by_counts.values, 0.02) # 2
OBJ = OBJ[(OBJ.obs.n_genes_by_counts < upper_lim) & (OBJ.obs.n_genes_by_counts > lower_lim)]
OBJ = OBJ[OBJ.obs.pct_counts_mt < 4] # 4

# Normalise every cell to 10,000 UMIs to account for difference in sequencing depth
sc.pp.normalize_total(OBJ, target_sum=1e4)
sc.pp.log1p(OBJ)

# Creates column 'highly_variable'
sc.pp.highly_variable_genes(OBJ, min_mean=0.0125, max_mean=3, min_disp=0.5)
OBJ = OBJ[:, OBJ.var.highly_variable]

sc.pp.regress_out(OBJ, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(OBJ, max_value=10)

OBJ.write_h5ad(
    'filtered_Pearsons1.h5ad',
    compression=hdf5plugin.FILTERS["zstd"]
)
