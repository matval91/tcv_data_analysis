#!/usr/bin/env python2.7


# routine to get values from a shot and plot them
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.interpolate as interp
shot=input('Shot number? ')
t=mds.Tree('tcv_shot', shot)
expdata={}
col=['k', 'r','b']

#plot 1
#ne autofit
ne_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:te')
ne = ne_s.data() #rho, time
ne_t = ne_s.getDimensionAt(0).data()
ne_rho = ne_s.getDimensionAt(1).data()
ne_err_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:te:error_bar')
ne_err = ne_err_s.data()
expdata[0] = {'x':ne_rho, \
              'y':ne, \
              't': ne_t, \
              'err': ne_err,\
              'ylab':r'n$_e$ $10^{20}\,(m^{-3})$'\
              }
print('Got te '+r'\tcv_shot::top.results.thomson.profiles.auto:te')


#ne_raw
ne_raw_s = t.getNode(r'\tcv_shot::top.results.thomson:te')
ne_raw = ne_raw_s.data()
ne_raw_t = ne_raw_s.getDimensionAt(0).data()
ne_raw_err_s = t.getNode(r'\tcv_shot::top.results.thomson:te:error_bar')
ne_raw_err = ne_raw_err_s.data()

polflux_th=t.getNode(r'\tcv_shot::top.results.thomson:psiscatvol')
polflux_time = polflux_th.getDimensionAt(0).data()
polflux_th = polflux_th.data() #R, time
polflux_n_th = np.copy(polflux_th)
rho_th = np.copy(polflux_th)
#finding rhopol
for i, el in enumerate(polflux_time):
    polflux_n_th[:, i] = (polflux_th[:,i]-max(polflux_th[:,i]))/(min(polflux_th[:,i])-max(polflux_th[:,i]))
    rho_th[:,i] = np.sqrt(polflux_n_th[:,i])
    
expdata[1] = {'x': rho_th, \
              'y': ne_raw, \
              't': ne_raw_t, \
              'err': ne_raw_err,\
              'ylab':r'n$_e$ $10^{20}\,(m^{-3})$'\
              }
print('Got te '+r'\tcv_shot::top.results.thomson:te')


# plot
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=10)
plt.rc('figure', facecolor='white')

ncol=4; nrow=4
ntot=ncol*nrow
timearr = np.linspace(0.4, min(expdata[0]['t'][-1], expdata[1]['t'][-1]), ntot)
f, ax = plt.subplots(nrow, ncol, sharex='all', sharey='all', figsize=[3.5*ncol,2.5*nrow])
f.suptitle('Shot #'+str(shot), fontsize=20)
for i in range(ncol):
    for j in range(nrow):
        ind=i*nrow+j

        #plot fit
        indfit = np.argmin(expdata[0]['t']-timearr[ind]<0.)
        x = expdata[0]['x']
        y = expdata[0]['y'][:, indfit]
        err = expdata[0]['err'][:, indfit]
        ax[i,j].errorbar(x,y,yerr=err,color='k', lw=2., label='fit '+str(expdata[0]['t'][indfit]))

        #plot raw
        indraw = np.argmin(expdata[1]['t']-timearr[ind]<0.)
        x = expdata[1]['x'][:, indraw]
        y = expdata[1]['y'][:, indraw]
        err = expdata[1]['err'][:, indfit]

        ax[i,j].errorbar(x,y,yerr=err, color='r', lw=2., label='raw '+str(expdata[1]['t'][indraw]))
        ax[i,j].legend(loc='best')
        ax[i,j].grid('on')
  
ax[0,0].set_xlim([0., 1])
for j in range(ncol):
    ax[nrow-1,j].set_xlabel(r'$\rho_{pol}$');
for i in range(nrow):
    ax[i,0].set_ylabel(r'$T_e$ (eV)');
plt.tight_layout()
plt.subplots_adjust(top=0.95)
plt.show()
