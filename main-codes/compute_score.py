

from anndata import AnnData
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import scdrs
import warnings
import scanpy as sc
import gc
from sklearn.externals import joblib


ff1= ['%s','%.5f','%.8f',"%.4e","%.4e",'%.5f','%.6f']
ff2=['%.5f']*1000
ff = ff1 + ff2

def df2csv(df,fname,myformats=[],sep='\t'):
  """
    # function is faster than to_csv
    # 7 times faster for numbers if formats are specified, 
    # 2 times faster for strings.
    # Note - be careful. It doesn't add quotes and doesn't check
    # for quotes or separators inside elements
    # We've seen output time going down from 45 min to 6 min 
    # on a simple numeric 4-col dataframe with 45 million rows.
  """
  if len(df.columns) <= 0:
    return
  Nd = len(df.columns)
  Nd_1 = Nd - 1
  formats = myformats[:] # take a copy to modify it
  Nf = len(formats)
  # make sure we have formats for all columns
  if Nf < Nd:
    for ii in range(Nf,Nd):
      coltype = df[df.columns[ii]].dtype
      ff = '%s'
      if coltype == np.int64:
        ff = '%d'
      elif coltype == np.float64:
        ff = '%f'
      formats.append(ff)
  fh=open(fname,'w')
  fh.write('\t'.join(df.columns) + '\n')
  for row in df.itertuples(index=False):
    ss = ''
    for ii in range(Nd):
      ss += formats[ii] % row[ii]
      if ii < Nd_1:
        ss += sep
    fh.write(ss+'\n')
  fh.close()




#OUT_FOLDER='$mydir/$myprojrct/output.scDRS'
#df_gs = scdrs.util.load_gs("$mydir/$myprojrct/zscore_merged2.gs")
#ls = pd.read_csv("$mydir/$myprojrct/zscore_merged2.gs", sep="\t", index_col=0)




OUT_FOLDER = sys.argv[0]
gsfile = argv[1]
humanatlas_h5ad = argv[2]



from sklearn.externals import joblib
df_gs = scdrs.util.load_gs(gsfile)
adata=sc.read_h5ad(humanatlas_h5ad)
df_cov = pd.DataFrame(adata.obs["nFeature_RNA"])
adata.var.index=adata.var["features"]
scdrs.preprocess(adata,cov=df_cov)
# joblib.dump(adata, '$mydir/adata_full.preprocessed.pkl') 

ls = pd.read_csv(gsfile, sep="\t", index_col=0)
for i in range(0, n):
    # adata = joblib.load('$mydir/adata_full.preprocessed.pkl') 
    trait=ls.index[i]
    gene_list = df_gs[trait][0]
    gene_weight = df_gs[trait][1]
    df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000, return_ctrl_norm_score=True)
    df_res.insert(0, "cell_id", df_res.index)
    df2csv(df_res,
           os.path.join(OUT_FOLDER, "%s.full_score" % trait),
           myformats=ff) 
    
    df2csv(df_res.iloc[:, 0:7],
           os.path.join(OUT_FOLDER, "%s.score" % trait),
           myformats=ff1) 
    del df_res
    del adata
    gc.collect()

# merge norm BPS table 

import os
import pandas as pd

os.chdir(OUT_FOLDER)
df_gs = pd.read_csv("list", header=None) 
# "list" is a txt file for all outputfile names

dict_df_stats = {
    trait: pd.read_csv(trait, sep="\t", index_col=0)
    for trait in df_gs.iloc[ : ,0]
    }
trait_list = list(dict_df_stats.keys())

df_norm = pd.concat(
    [dict_df_stats[trait]["norm_score"] for trait in trait_list], axis=1
)
df_norm.columns = trait_list

df_norm.to_csv("1_df_norm_score.tsv", sep="\t")