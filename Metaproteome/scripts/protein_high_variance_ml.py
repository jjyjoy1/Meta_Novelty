import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import joblib

# Inputs
abundance_file = snakemake.input["abundance"]
metadata_file = snakemake.input["metadata"]
clusters_out = snakemake.output["clusters"]
model_out = snakemake.output["model"]

# Load data
df = pd.read_csv(abundance_file, index_col=0)
meta = pd.read_csv(metadata_file, index_col=0)

# Standardize
X = StandardScaler().fit_transform(df.T)

# High variance proteins
variances = df.var(axis=1)
top_proteins = variances.sort_values(ascending=False).head(200).index
X_highvar = df.loc[top_proteins].T

# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_highvar)

# Clustering
kmeans = KMeans(n_clusters=3, random_state=42).fit(X_pca)
clusters = pd.DataFrame(kmeans.labels_, index=X_highvar.index, columns=["Cluster"])
clusters.to_csv(clusters_out)

# ML: metadata prediction
y = meta.loc[X_highvar.index]["Condition"]  # Modify to your metadata column
clf = RandomForestClassifier(random_state=42)
clf.fit(X_highvar, y)
joblib.dump(clf, model_out)


