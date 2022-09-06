import numpy as np
import pandas as pd

with open('id.txt') as f:
    lines = f.readlines()

id = int(lines[0])
df = pd.read_csv('input_grid.csv', header=None)

with open('cont.txt') as f:
    lines = f.readlines()

cont = int(lines[0])

with open('tinit.txt') as f:
    lines = f.readlines()

tinit = float(lines[0])

with open('input.txt', 'w') as f:
        f.write(str(df.iloc[id,0]))
        f.write('\n')
        f.write(str(df.iloc[id,1]))
        f.write('\n')
        f.write(str(df.iloc[id,2]))
        f.write('\n')
        f.write(str(cont))
        f.write('\n')
        f.write(str(tinit))
        f.write('\n')