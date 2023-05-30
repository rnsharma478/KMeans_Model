import pandas as pd
import numpy as np 

dataset = pd.read_csv("/scf-data/kopal/parameter_map.csv", index_col = 0)
lst = ['a','b','c','d','e','f','g','h','i','j','k','l','ma','n','o','p','q','r','s','t','u','v','w','x','y','z','aa','ab','ac','ad','ae']

para_dict = {}
for z in lst:
  zero_arr = np.zeros((20,975))
  para_dict[z] = pd.DataFrame(zero_arr)

  for i in range(20):
    for j in range(31):
      for k in range(975):
        para_dict[z].iloc[i][k] = (list((dataset.iloc[j][i][1:-1]).split(",")))[k]

for z in lst:
  para_dict[z] = para_dict[z].T

import pandas as pd
from sklearn.cluster import KMeans

#X, _ = make_blobs(n_samples=10, centers=2, n_features=975)

#df = pd.DataFrame(X, columns=['Feat_1', 'Feat_2', 'Feat_3', 'Feat_4'])

y_dict = {}

kmeans = KMeans(n_clusters=2)

for z in lst:
  y = kmeans.fit_predict(para_dict[z])
  y_dict[z] = y

clusters = pd.DataFrame.from_dict(y_dict)
clusters.to_csv("clusters.csv")
