import numpy as np
import pandas as pd

ins = []

for i in np.arange(3.00,8.50,0.50):
    for j in np.arange(3.0,15.5,1.0):
        #for k in np.arange(0,2,1):
        ins.append([round(i,2),round(j,2),1])

# for i in np.arange(2.0,4.5,0.5):
#     for j in np.arange(0.5,10.5,0.5):
#         #for k in np.arange(0,2,1):
#         ins.append([round(i,2),round(j,2),1])


data = {'size': [i[0] for i in ins], 'massloss':[i[1] for i in ins], 'geom':[i[2] for i in ins]}

df = pd.DataFrame(data=data)

df.to_csv('input_grid_Al2O3_sept23.csv', header=False, index=False)