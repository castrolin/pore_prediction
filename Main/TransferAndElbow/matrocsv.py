import numpy as np 
import h5py
import xlsxwriter
import os
#import pandas as pd

def matfile_to_xlsx(mat_name,number):

    with h5py.File(mat_name,'r') as f:

        '''
        for model in f.keys():
            print(model)
            for i in f['DataBase']:
                print(f['DataBase'][i].shape[0])
        '''

        #length array
        Length = []
        for index in range(0,f['DataBase']['Length'].shape[0]):
            Length.append(''.join([str(i) for i in f[f['DataBase']['Length'][index][0]].value]))
    
        #Width array
        Width = []
        for index in range(0,f['DataBase']['Width'].shape[0]):
            Width.append(''.join([str(i) for i in f[f['DataBase']['Width'][index][0]].value]))
        #Ratio array
        Ratio = []
        for index in range(0,f['DataBase']['Ratio'].shape[0]):
            Ratio.append(''.join([str(i) for i in f[f['DataBase']['Ratio'][index][0]].value]))
        #Angle array
        Angle = []
        for index in range(0,f['DataBase']['Angle'].shape[0]):
            Angle.append(''.join([str(i) for i in f[f['DataBase']['Angle'][index][0]].value]))
        # write into xlsx file
        workbook = xlsxwriter.Workbook('Meltpool'+str(number)+'.xlsx')
        worksheet = workbook.add_worksheet() #worksheet.write(row,column, value)
        for index in range(0,len(Length)):
            worksheet.write(index,0,float(Length[index][1:-1]))
            worksheet.write(index,1,float(Width[index][1:-1]))
            worksheet.write(index,2,float(Ratio[index][1:-1]))
            worksheet.write(index,3,float(Angle[index][1:-1]))
    workbook.close()
    

## read multi matfiles with name
#Path = r"X:\Castro\pore_prediction\database0619"
#Path = r"X:\Castro\pore_prediction\exp0715
Path = r"X:\Castro\DataBase0715"

#for file in os.listdir():
    #if file.endswith('.mat'):
files = os.listdir(Path)
for file in files:
    if file.endswith('.mat'):
        print(file)
        for i in range(1,len(file)):
            if file[i] == ".":
                number = []
                for n in range(0,i):
                    number.append(file[n])

                print(number)

                if len(number) == 1:
                    sum = 0
                    for ii in range(0,len(number)):
                        sum = sum + int(number[ii])
                elif len(number) == 2:
                    sum = 0
                    for ii in range(0,len(number)):
                        sum += int(number[ii])*(10**(len(number)-1-ii))
                elif len(number) == 3:
                    sum = 0
                    for ii in range(0,len(number)):
                        sum += int(number[ii])*(10**(len(number)-1-ii))
                print(sum)
                matfile_to_xlsx("X:\Castro\DataBase0715\\"+file,sum)



#print(files)








'''

How to extract the data from mat file and release the dataset!!
1. Build the code for python SOM algortithm based on csv file 
2. .mat to .csv and do the same things
3. matlab is too slow and non-sense!!
4. print(Length[139][1:-1]) #[1:-1] remove the string bracket

Solve the problem ^^

'''



