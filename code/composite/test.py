# %%
import os
import pandas as pd
import sccomposite
from sccomposite import RNA_modality, ADT_modality, ATAC_modality, Multiomics

# %%
data_dir = '../../data-raw/composite'
#out_dir = './rna'
out_dir = './rna_adt'
os.makedirs(out_dir, exist_ok=True)
rna_mtx_path = os.path.join(data_dir, 'RNA.mtx')
adt_mtx_path = os.path.join(data_dir, 'ADT.mtx')

# %%
#multiplet_classification, consistency = RNA_modality.composite_rna(rna_mtx_path)
#score = consistency

multiplet_classification, multiplet_probability = Multiomics.composite_multiomics(RNA=rna_mtx_path, ADT=adt_mtx_path)
score = multiplet_probability

# %%
#multiplet_classification = ["singlet", "doublet", "singlet"]
#score = [0.5, 0.5, 0.5]
data_dict = {'composite_classification': multiplet_classification,
             'score': score}
df = pd.DataFrame(data_dict)
print(df.head())
df.index.name = 'index'
df.reset_index(inplace=True)
df.to_csv(os.path.join(out_dir, 'composite_classification.tsv'), 
          index=False, 
          sep=str("\t"))
