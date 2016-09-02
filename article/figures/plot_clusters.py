
"""
Produce plots for the GCs (HRD, memebership, and chemistry). This code may be clunky.
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from astropy.table import join
from astropy.table import Table



#rave_cannon_dr1 = Table.read('/Users/khawkins/Desktop/RAVE_cannon/rave-tgas-v1.fits')
#RAVEDR4_GC = Table.read('/Users/khawkins/Desktop/RAVE_cannon/RAVEDR4_GC.fits')
#data_table = join(rave_cannon_dr1, RAVEDR4_GC, keys=("Name",))
#RAVEDR4 = Table.read('/Users/khawkins/Desktop/RAVE_cannon/RAVE-DR4.fits')
#DR4_cannon = join(rave_cannon_dr1, RAVEDR4, keys=("Name",))


try:
    data_table

except NameError:
    from rave_io import get_cannon_dr1 

    RAVEDR4_GC = Table.read('/data/gaia-eso/kh536/RAVEDR4_GC.fits')
    keep = ("Name", "NGC", "RAVEID")
    for column in RAVEDR4_GC.dtype.names:
        if column not in keep:
            del RAVEDR4_GC[column]

    data_table = join(get_cannon_dr1(), RAVEDR4_GC, keys=("Name",))

else:
    print("Warning: Using pre-loaded data!")

DR4_cannon = data_table


#---define unique clusters
clusters = np.unique(data_table['NGC'])
#-------Generate HRD for each cluster----
limits = [[3400,7000], [0.0,4.99]]
K = len(clusters)
factor = 3.5
lbdim = 0.2 * factor
trdim = 0.1 * factor
whspace = 0.05
yspace = factor
xspace = factor * K + factor * (K - 1) * whspace + lbdim * (K - 1)
xdim = lbdim + xspace + trdim
ydim = lbdim + yspace + trdim

fig, axes = plt.subplots(K,2, figsize=(xdim, ydim), sharex=True,sharey=True)
fig.subplots_adjust(left=lbdim/xdim, bottom=lbdim/ydim, 
  right=(xspace + lbdim)/xdim, top=(yspace + lbdim)/ydim, 
  wspace=whspace, hspace=whspace)


for i in np.arange(K):
  cind = np.where(data_table['NGC'] == clusters[i])[0] #find the cluster stars
  #plot the cannon
  ax = axes[i][0]
  vmin, vmax = -2, 0.5
  scat = ax.scatter(data_table['TEFF'][cind], data_table['LOGG'][cind], c=data_table['FE_H'][cind], vmin=vmin, vmax=vmax)
  ax.set_xlim(limits[0])
  ax.set_ylim(limits[1])
  ax.xaxis.set_major_locator(MaxNLocator(6))
  ax.yaxis.set_major_locator(MaxNLocator(6))
  #ax.set_axis_bgcolor("#CCCCCC")
  ax.invert_xaxis();ax.invert_yaxis()
  ax.text(6800,1.0,'NGC%i'%clusters[i])
  ax.set_ylabel(r"$\log{g}$",)

  if i==0:
  	ax.set_title('The Cannon')
  if i == K-1: 
    ax.set_xlabel(r"$T_{\rm eff}$")

  ax = axes[i][1]
  #plot the RAVE DR4
  ax.scatter(data_table['TeffK'][cind], data_table['loggK'][cind], c=data_table['__M_H_K'][cind], vmin=vmin, vmax=vmax)
  ax.set_xlim(limits[0])
  ax.set_ylim(limits[1])
  ax.xaxis.set_major_locator(MaxNLocator(6))
  ax.yaxis.set_major_locator(MaxNLocator(6))
  #ax.set_axis_bgcolor("#CCCCCC")
  ax.invert_xaxis();ax.invert_yaxis()
  if i==0:
  	ax.set_title('RAVE DR4')
  if i == K-1: 
    ax.set_xlabel(r"$T_{\rm eff}$") 
    ax.get_yaxis().set_visible(False)

plt.colorbar(scat)

def generate_cluster_abundance_plot(cluster=clusters[0]):
  all_columns = [ 
  ("C_H", "N_H"), 
  ("C_H", "NA_H"),
  ]
  cind = np.where(data_table['NGC'] == cluster)[0] #find the cluster stars
  RAlim=5
  feild_ind = np.where((DR4_cannon['RAJ2000'] > data_table['RAJ2000'][cind[0]] - RAlim)&(DR4_cannon['RAJ2000'] < data_table['RAJ2000'][cind[0]] + RAlim) &\
    (DR4_cannon['DEJ2000'] > data_table['DEJ2000'][cind[0]] - RAlim)&(DR4_cannon['DEJ2000'] < data_table['DEJ2000'][cind[0]] + RAlim))[0]

  plt.figure()
  plt.plot(DR4_cannon['RAJ2000'][feild_ind],DR4_cannon['DEJ2000'][feild_ind],'.',color='gray')
  plt.plot(data_table['RAJ2000'][cind],data_table['DEJ2000'][cind],'bo',ms=7)


  fig, axes = plt.subplots(len(all_columns), 1)
  for i, (ax, columns, ) in enumerate(zip(axes, all_columns,)):
    x = data_table[columns[0]][cind]-data_table['FE_H'][cind]
    y = data_table[columns[1]][cind]-data_table['FE_H'][cind]
    xfield = DR4_cannon[columns[0]][feild_ind]
    yfield = DR4_cannon[columns[1]][feild_ind]
    ax.plot(xfield,yfield,'.',color='gray',label='Field',ms=2)
    ax.plot(x,y,'bo', label='NGC%i'%cluster,ms=7)
    ax.set_ylim([-1,1]);ax.set_xlim([-1,1])
    #ax.legend()
    ax.set_xlabel('[%s/Fe]'%columns[0].split('_')[0])
    ax.set_ylabel('[%s/Fe]'%columns[1].split('_')[0])
    if i == 0:
      ax.legend(fontsize='small')
  plt.tight_layout()


def generate_cluster_memebership(cluster=clusters[0]):
  all_columns = [ 
  ("HRV", "FE_H"), 
  ("HRV", "__M_H_K"),
  ]
  cind = np.where(data_table['NGC'] == cluster)[0] #find the cluster stars
  RAlim=5
  feild_ind = np.where((DR4_cannon['RAJ2000'] > data_table['RAJ2000'][cind[0]] - RAlim)&(DR4_cannon['RAJ2000'] < data_table['RAJ2000'][cind[0]] + RAlim) &\
    (DR4_cannon['DEJ2000'] > data_table['DEJ2000'][cind[0]] - RAlim)&(DR4_cannon['DEJ2000'] < data_table['DEJ2000'][cind[0]] + RAlim))[0]

  x = data_table['HRV'][cind]
  y = data_table['FE_H'][cind]
  xfield = DR4_cannon['HRV'][feild_ind]
  yfield = DR4_cannon['FE_H'][feild_ind]
  plt.figure(figsize=(xdim, ydim))
  ax1 = plt.subplot(1,2,1)
  ax1.plot(xfield,yfield,'.',color='gray',label='Field',ms=2)
  ax1.plot(x,y,'bo',label='NGC%i'%cluster,ms=7)
  ax1.set_title('The Cannon')
  ax1.set_ylabel('[Fe/H]')
  ax1.set_xlabel('HRV (km/s)')
  ax1.set_ylim([-2.0,0.6]) ; ax1.set_xlim([-400,400])
  ax1.legend(fontsize='small')



  ax2 = plt.subplot(1,2,2,sharex=ax1,sharey=ax1)
  y = data_table['__M_H_K'][cind]
  xfield = DR4_cannon['HRV'][feild_ind]
  yfield = DR4_cannon['__M_H_K'][feild_ind]
  ax2.plot(xfield,yfield,'.',color='gray',label='Field',ms=2)
  ax2.plot(x,y,'bo',label='NGC%i'%cluster,ms=7)
  ax2.set_title('RAVE DR4')
  ax2.set_ylabel('[Fe/H]')
  ax2.set_xlabel('HRV (km/s)')
  ax2.set_ylim([-2.0,0.6]) ; ax2.set_xlim([-400,400])     
  ax2.get_yaxis().set_visible(False)

  plt.tight_layout()



for i in np.arange(len(clusters)):
  #generate_cluster_abundance_plot(cluster=clusters[i])
  generate_cluster_memebership(cluster=clusters[i])






plt.show()




