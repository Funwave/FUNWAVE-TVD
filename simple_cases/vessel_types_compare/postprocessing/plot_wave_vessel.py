### Plot Flat Vessel Waves ###

# import necessary modules
import os
import numpy as np               
import matplotlib.pyplot as plt


# write your OWN PC folder path for fdir
# Remember that we use for Mac & Linux machines '/', while on windows '\'
HOME = os.environ['HOME']
fdir = os.path.join(HOME,'FUNWAVE-TVD','simple_cases','vessel_flat_bottom','work', 'output')
fileName = os.path.join(fdir,'eta_00001')
eta = np.loadtxt(fileName)

# define plot location
n,m = np.shape(eta)
dx = 1.0
dy = 1.0

x = np.asarray([float(xa)*dx for xa in range(m)])
y = np.asarray([float(ya)*dy for ya in range(n)])

nfile = [20, 40]    # range of eta files you want to plot
sec = ['20','40']  # time  you want to plot

# figure size option 
wid = 8    # width
length = 5 # length

# Plot figure
fig = plt.figure(figsize=(wid,length),dpi=200)

for num in range(len(nfile)):
    fnum= '%.5d' % nfile[num]
    etaFile = os.path.join(fdir,'eta_'+fnum)
    eta = np.loadtxt(etaFile)

    ax = fig.add_subplot(len(nfile),1,num+1)
    fig.subplots_adjust(hspace=.45)
    plt.pcolor(x, y, eta,cmap='coolwarm')
    
    title = 'Time = '+sec[num]+ ' sec'
    plt.title(title)
    plt.axis('tight')

    plt.ylabel('Y (m)')
    plt.xlabel('X (m)')
    cbar=plt.colorbar()
    cbar.set_label(r'$\eta$'+' (m)', rotation=90)
    plt.clim(-1.5, 1.5)

fig.savefig('eta_flat_vessel.png', dpi=fig.dpi)
