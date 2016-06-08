import numpy as np
import matplotlib, datetime
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation


import pandas as pd
import sys,os,re,multiprocessing,netCDF4
from netCDF4 import Dataset

#netcdf file path
ncfile = sys.argv[1]

netCDF_data = Dataset(ncfile, mode='r')
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


def colourmap(mn=0,mx=1):
     #returns function that when used between range mn and mx give the rgb equivalent of the colorscheme selected
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    jet = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=mn, vmax=mx)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    return scalarMap.to_rgba




def update_plt(num):
    timestep = runs[num]
    M=s_file.loc[timestep,'M']       
    shead = s_file.columns
    if  num%10 == 0 : print num
    cm(num)
    plt.cla()
    plt.title(datetime.datetime.fromtimestamp(timestep*10*60).strftime('%d %H:%M'))
    #plt.xlim(0, 1)
    plt.ylim(-200, 100)
    plt.xlabel('index')
    plt.ylabel('conc/M')
    
    line = [plt.plot(s_file.ix[timestep].map(lambda x: np.log10(x/s_file.loc[timestep,'M'])) ,color=cm(num), label=num)]#,'r'

    return line,



# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


for group in netCDF_data.groups:
    print '\n\nReading', group

    global s_file,cm,runs
    s_file = pd.DataFrame(netCDF_data.groups[group].variables['Spec'][:])
    s_file.columns = str(netCDF_data.groups[group].variables['Spec'].head).split(',')
    
    runs = s_file.index# [0] #### LAST RUN NEEDS TO BE ONE WITH GREATEST NO REACTIONS! ### 144,288,431,    
    cm = colourmap(min(runs),max(runs))
    fig1 = plt.figure()

    line_ani = animation.FuncAnimation(fig1, update_plt, len(runs), fargs=(), interval=50, blit=True)
    
    line_ani.save('%s.mp4'%group, writer=writer)








