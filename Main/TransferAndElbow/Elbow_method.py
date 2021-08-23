import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import pandas as pd
from numba import jit
import os
import xlsxwriter  

#importing the dataset

# Using the elbow method to find  the optimal number of clusters

clusternumber = []

@jit
def elbow_method(Maimum_cluster,dataset):
  wcss =[]
  for i in range(1, Maimum_cluster):
    kmeans=KMeans(n_clusters=i, init='k-means++', max_iter= 300, n_init= 10, random_state= 0)
    kmeans.fit(dataset)
    wcss.append(kmeans.inertia_)
  print('.....EXCUTE....')
  return wcss

@jit
def number_of_cluster(wcss, Max_percentage,Min_percentage):
  print('.....EXCUTE....')
  Percentage_of_wcss = []
  Percentage_of_variance = []
  maximum_wcss = max(wcss)
  for i in range(0,len(wcss)):
    Percentage_of_wcss.append((wcss[i]/maximum_wcss)*100)
    Percentage_of_variance.append(1-Percentage_of_wcss[i])

  n = 0
  for numberofcluster in range(0,len(Percentage_of_wcss)):
    if Percentage_of_wcss[numberofcluster] <= Max_percentage and Percentage_of_wcss[numberofcluster] > Min_percentage :
      n +=1
      if n == 1:
        print(numberofcluster)
        return numberofcluster
  #plt.plot(range(1,len(Percentage_of_wcss)+1),Percentage_of_wcss)
  #plt.show()


#...................................
Path = r"X:\Castro\Data_transform"
files = os.listdir(Path)
n = 0
# create the excel file to save the cluster number of each excel
workbook = xlsxwriter.Workbook('Number_of_cluster.xlsx')
worksheet = workbook.add_worksheet()
filenumber = []
for file in files:
  if file.endswith('.xlsx'):
    #print(file)
    for i in range(7,len(file)): # M e l t p o o l is 7 array
      if file[i] == ".":
        number = []
        for n in range(8,i):
          number.append(file[n])
        #print(number)

        if len(number) == 1:
          sum =0
          for ii in range(0,len(number)):
            sum += int(number[ii])
        elif len(number) == 2:
          sum = 0
          for ii in range(0,len(number)):
            sum += int(number[ii])*(10**(len(number)-1-ii))
        elif len(number) == 3:
          sum = 0
          for ii in range(0,len(number)):
            sum += int(number[ii])*(10**(len(number)-1-ii))
        print(sum)
        filenumber.append(sum)
    # Excute the elbow method and data_processingg

    dataset=pd.read_excel(file,engine='openpyxl')
    print('loading file:',file)
    X=dataset.iloc[: , :].values
    x1=dataset.iloc[:,0]
    y1=dataset.iloc[: , 2:3]
    y1=y1.values.reshape(-1, 1)
    #print(pd.DataFrame(dataset))
    #print(X)

    #Excute the result
    wcss = elbow_method(50,X)
    numberofCluster = number_of_cluster(wcss,10,9)
    print('number of Cluster:',numberofCluster)
    clusternumber.append(numberofCluster)
    del X
    del wcss
    del numberofCluster
    n+=1

# Check the list: Did it has a None Value
clusternumber_mean = 0
n = 0
for index in range(0, len(clusternumber)):
  if clusternumber[index] != None:
    clusternumber_mean += clusternumber[index]
    n += 1
clusternumber_mean = clusternumber_mean/n
for index in range(0,len(clusternumber)):
  if clusternumber[index] == None:
    clusternumber[index] = int(clusternumber_mean)



# Save to the excel file    
for num in range(0,len(clusternumber)):
  worksheet.write(num,1,clusternumber[num])
  worksheet.write(num,0,filenumber[num])
workbook.close()

print(clusternumber)
plt.plot(clusternumber) 
plt.show()    
    




#...................................

'''  
wcss = elbow_method(30,X)
print(wcss)
maximum_wcss = max(wcss)
Percentage_of_wcss = (wcss/maximum_wcss)*100
print(Percentage_of_wcss)
plt.plot(range(1, 30),Percentage_of_wcss)
plt.title('The Elbow Method')
plt.xlabel('Number of clusters K')
plt.ylabel('Average Within-Cluster distance to Centroid (WCSS) %')  
plt.show()


''' 