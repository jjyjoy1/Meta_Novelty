import pandas as pd
from scipy.stats import zscore

input_file = snakemake.input[0]
output_file = snakemake.output[0]
sample_name = snakemake.wildcards.sample

df = pd.read_csv(input_file)
df.set_index("Protein IDs", inplace=True)

z_scores = zscore(df, nan_policy='omit', axis=1)
sample_scores = pd.Series(z_scores[:, df.columns.get_loc(sample_name)], index=df.index)

outliers = sample_scores[abs(sample_scores) > 2].sort_values(ascending=False)
outliers.to_csv(output_file, header=["Z-score"])



