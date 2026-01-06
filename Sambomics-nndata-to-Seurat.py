## This script converts h5ad anndata to seurat-readable files.
## The script has been copied from sanbomics_scripts (mousepixels - github).

import scanpy as sc
from scipy import io
!mkdir matrix_files
!mkdir 
adata = sc.read_h5ad("zf_atlas_24hpf_v1_release.h5ad") 
adata

with open('matrix_files/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')
        

with open('matrix_files/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')
        
io.mmwrite('matrix_files/matrix', adata.X.T)


!ls matrix_files/
!gzip matrix_files/*
