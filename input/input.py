import numpy as np
import pandas as pd

df = pd.read_csv('input_grid.csv', header=None)

with open('run.in') as f:
    lines = f.readlines()

id    = int(lines[0].strip('\n'))
cont  = int(lines[1].strip('\n'))
tinit = float(lines[2].strip('\n'))
tend  = float(lines[3].strip('\n'))
composition = lines[4].strip('\n')
sdist = int(lines[5].strip('\n'))
dist_type = int(lines[6].strip('\n'))
mu    = float(lines[7].strip('\n'))
sigma = float(lines[8].strip('\n'))
mdot = float(lines[9].strip('\n'))
geom=int(lines[10].strip('\n'))

if (dist_type==1):
    df_dist = pd.read_csv('distributions.in', header=None, skiprows=1)
if (dist_type==2):
    df_dist = pd.read_csv('n_pow_dist.dat', header=None)

with open('input.txt', 'w') as f:
        if (dist_type>=1):
            f.write(str(1.0))
            f.write('\n')
            f.write(str(mdot))
            f.write('\n')
            f.write(str(geom))
        else:
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
        f.write(str(tend))
        f.write('\n')
        f.write(str(sdist))
        f.write('\n')
        f.write(str(dist_type))
        f.write('\n')
        if (dist_type==1):
            f.write(str(df_dist.iloc[id,0]))
            f.write('\n')
            f.write(str(df_dist.iloc[id,1]))
            f.write('\n')

        elif (dist_type==2):
            f.write(str(df_dist.iloc[id,0]))
            f.write('\n')
            f.write(str(sigma))
            f.write('\n')
        else:
            f.write(str(mu))
            f.write('\n')
            f.write(str(sigma))
            f.write('\n')
        
        f.write(composition)
        f.write('\n')
        
