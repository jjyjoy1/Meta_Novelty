import pandas as pd
from sklearn.impute import KNNImputer
from sklearn.decomposition import TruncatedSVD

abundance_file = snakemake.input["abundance"]
output_file = snakemake.output["novel"]

# Load data
df = pd.read_csv(abundance_file, index_col=0)
df.replace(0, pd.NA, inplace=True)

# Focus on low abundance: proteins in bottom 25% mean intensity
low_mean_proteins = df.mean(axis=1).nsmallest(int(len(df)*0.25)).index
low_df = df.loc[low_mean_proteins]

# Imputation with matrix completion
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(low_df.T)
imputed_df = pd.DataFrame(imputed_data.T, index=low_df.index, columns=low_df.columns)

# Detect novel proteins: high variance across imputed but low mean original
diff = imputed_df.mean(axis=1) - low_df.mean(axis=1)
novel_candidates = diff.sort_values(ascending=False).head(20)
novel_candidates.to_csv(output_file, header=["novelty_score"])


