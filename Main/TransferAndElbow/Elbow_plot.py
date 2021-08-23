import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import pandas as pd
from numba import jit
import os
import xlsxwriter  
file = r"X:\Castro\Data_transform\Meltpool63.xlsx"
dataset=pd.read_excel(file,engine='openpyxl')
print('loading file:',file)
X=dataset.iloc[: , :].values
x1=dataset.iloc[:,0]
y1=dataset.iloc[: , 2:3]
y1=y1.values.reshape(-1, 1)

wcss = []
for i in range(1,50):
    kmeans=KMeans(n_clusters=i, init='k-means++', max_iter= 300, n_init= 10, random_state= 0)
    kmeans.fit(dataset)
    wcss.append(kmeans.inertia_)

Percentage_of_wcss = []
Percentage_of_variance = []
maximum_wcss = max(wcss)
for i in range(0,len(wcss)):
    Percentage_of_wcss.append((wcss[i]/maximum_wcss)*100)
    Percentage_of_variance.append(100-Percentage_of_wcss[i])

n = 0
for numberofcluster in range(0,len(Percentage_of_wcss)):
    if Percentage_of_wcss[numberofcluster] <= 10 and Percentage_of_wcss[numberofcluster] > 9 :
      n +=1
      if n == 1:
        print(numberofcluster)

plt.plot(range(1,len(Percentage_of_wcss)+1),Percentage_of_variance)
plt.show()