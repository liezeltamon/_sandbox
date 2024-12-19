# %%
import os
import pandas as pd
import sccomposite
from sccomposite import RNA_modality, ADT_modality, ATAC_modality, Multiomics

# %%
data_dir = '../../data-raw/composite'
out_dir = os.makedirs('./rna', exist_ok=True)
rna_mtx_path = os.path.join(data_dir, 'RNA.mtx')
adt_mtx_path = os.path.join(data_dir, 'ADT.mtx')

# %%
multiplet_classification, consistency = RNA_modality.composite_rna(rna_mtx_path)
reliability = consistency

#multiplet_classification, reliability = Multiomics.composite_multiomics(RNA = rna_mtx_path, ADT =  adt_mtx_path)
#reliability = reliability

# %%
data = {'multiplet_classification': multiplet_classification}
data_file = pd.DataFrame(data)
print(data_file.head())

# Save reliability to CSV
probability_df = pd.DataFrame({'reliability': reliability})
print(probability_df.head())

#probability_df.index.name = 'index'
#probability_df.reset_index(inplace=True)
#probability_df.to_csv(os.path.join(out_dir, "reliability.csv"), index=False)

data_file.index.name = 'index'
data_file.reset_index(inplace=True)
data_file.to_csv(out_dir, "multiplet_prediction.csv", index=False, sep = str(','))
