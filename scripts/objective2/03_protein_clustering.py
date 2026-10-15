from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


INPUT = Path(
    "data/processed/annotation/alternaria/protein_features.csv"
)
OUTPUT_DIR = Path(
    "data/processed/annotation/alternaria"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

# Select numerical feature columns; exclude the protein ID.
feature_columns = [
    col for col in df.columns
    if col != "protein_id"
]

X = df[feature_columns].replace(
    [float("inf"), float("-inf")], float("nan")
)

# Exclude rows with invalid feature values.
valid = X.notna().all(axis=1)
df = df.loc[valid].reset_index(drop=True)
X = X.loc[valid].reset_index(drop=True)

if len(df) < 4:
    raise ValueError("At least 4 valid proteins are needed.")

# Standardize the numerical features.
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Exploratory clustering.
model = KMeans(
    n_clusters=4,
    random_state=42,
    n_init=10
)

df["cluster"] = model.fit_predict(X_scaled)

# Two-dimensional PCA projection for plotting.
pca = PCA(n_components=2)
coordinates = pca.fit_transform(X_scaled)

df["PCA1"] = coordinates[:, 0]
df["PCA2"] = coordinates[:, 1]

# Save protein-level cluster assignments.
df.to_csv(
    OUTPUT_DIR / "protein_clusters.csv",
    index=False
)

# Save a cluster summary.
summary = df.groupby("cluster")[feature_columns].mean()
summary["protein_count"] = df.groupby("cluster").size()

summary.to_csv(
    OUTPUT_DIR / "cluster_summary.csv"
)

# Create a visualization.
plt.figure(figsize=(10, 7))

for cluster_id in sorted(df["cluster"].unique()):
    subset = df[df["cluster"] == cluster_id]

    plt.scatter(
        subset["PCA1"],
        subset["PCA2"],
        s=10,
        alpha=0.6,
        label=f"Cluster {cluster_id}"
    )

plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.title("Exploratory Clustering of A. alternata Proteins")
plt.legend()
plt.tight_layout()

plt.savefig(
    OUTPUT_DIR / "protein_clusters_pca.png",
    dpi=200
)
plt.close()

print("Clustering completed.")
print(f"Proteins analyzed: {len(df)}")
print(f"Cluster counts:\n{df['cluster'].value_counts().sort_index()}")
print(f"Explained variance by PCA: {pca.explained_variance_ratio_}")