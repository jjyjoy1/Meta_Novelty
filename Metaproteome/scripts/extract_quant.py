import pandas as pd

input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Load proteinGroups.txt
df = pd.read_csv(input_file, sep='\t')

# Extract columns with LFQ intensities
lfq_columns = [col for col in df.columns if col.startswith("LFQ intensity ")]
protein_ids = df["Protein IDs"]

# Form the abundance matrix
lfq_data = df[lfq_columns]
lfq_data.columns = [col.replace("LFQ intensity ", "") for col in lfq_columns]
lfq_data.insert(0, "Protein IDs", protein_ids)

# Save output
lfq_data.to_csv(output_file, index=False)


