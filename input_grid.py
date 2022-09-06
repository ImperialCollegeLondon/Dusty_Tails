import numpy as np
import pandas as pd

ins = []

for i in np.arange(1.0,5.5,0.5):
    for j in np.arange(1.0,10.5,1.0):
        for k in np.arange(0,2,1):
            ins.append([i,j,k])


data = {'size': [i[0] for i in ins], 'massloss':[i[1] for i in ins], 'geom':[i[2] for i in ins]}

df = pd.DataFrame(data=data)

df.to_csv('input_grid.csv', header=False, index=False)