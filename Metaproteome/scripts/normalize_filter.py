import pandas as pd
import numpy as np

input_file = snakemake.input[0]
output_file = snakemake.output[0]

df = pd.read_csv(input_file)

# Replace 0s with NaNs
df.replace(0, np.nan, inplace=True)

# Drop proteins with more than 50% missing values
df = df.dropna(thresh=int(0.5 * df.shape[1]))

# Log2 transform
df.iloc[:, 1:] = np.log2(df.iloc[:, 1:])

# Quantile normalization
rank_mean = df.iloc[:, 1:].stack().groupby(df.iloc[:, 1:].rank(method='first').stack().astype(int)).mean()
df.iloc[:, 1:] = df.iloc[:, 1:].rank(method='min').stack().astype(int).map(rank_mean).unstack()

# Save
df.to_csv(output_file, index=False)

