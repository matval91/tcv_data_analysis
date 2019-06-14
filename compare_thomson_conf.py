# routine to get values from a shot and plot them
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.interpolate as interp
import eqtools

#shot=input('Shot number? ')
shot=58823
t=mds.Tree('tcv_shot', shot)
expdata={}
expprof={}
col=['k', 'r','b']
eq=eqtools.TCVLIUQETree(shot)

################
#THOMSON DATA
################

#ne
ne_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:ne')
rho_thomson_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:rho') #normalized poloidal flux
rho_thomson = rho_thomson_s.data()**0.5
ne = ne_s.data()
ne_t = ne_s.getDimensionAt(0).data()
ne0 = ne[0,:]
ne5 = ne[20,:]
expdata[0] = {'nq':2, 'x':[ne_t, ne_t], \
              'y':[ne0*1e-20, ne5*1e-20], 'ylab':r'n$_e$ $10^{20}\,(m^{-3})$',\
              'lab':[r'$\rho_{POL}=0$', r'$\rho_{POL}=0.5$']}
expprof[0] = { 'x':rho_thomson, 't':ne_t, 'y':ne*1e-20, 'ylab':r'n$_e$ $10^{20}\,(m^{-3})$',\
               'lab':'thomson'}

print('Got ne from thomson')

#plot 7
#te
te_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:te')
te = te_s.data()
te_t = te_s.getDimensionAt(0).data()
te0 = te[0,:]
te5 = te[20,:]
expdata[1] = {'nq':2, 'x':[te_t, te_t], \
              'y':[te0, te5], 'ylab':r'T$_e$ (keV)',\
              'lab':[r'$\rho_{POL}=0$', r'$\rho_{POL}=0.5$']}
expprof[1] = { 'x':rho_thomson, 't':te_t, 'y':te*1e-3, 'ylab':r'T$_e$ (keV)$',\
               'lab':'thomson'}

print('Got te from thomson')



########################
# CONF DATA
########################
ne_s=t.getNode(r'\tcv_shot::top.results.conf:ne')
rho_s=t.getNode(r'\tcv_shot::top.results.conf:rhotor')

rho_conf_ne = rho_s.data()
ne = ne_s.data()
ne_t = ne_s.getDimensionAt(1).data()
ne0 = ne[:,0]
ne5 = ne[:,20]
expdata[2] = {'nq':2, 'x':[ne_t, ne_t], \
              'y':[ne0*1e-20, ne5*1e-20], 'ylab':r'n$_e$ $10^{20}\,(m^{-3})$',\
              'lab':[r'$\rho_{TOR}=0$', r'$\rho_{TOR}=0.5$']}
expprof[2] = { 'x':rho_conf_ne, 't':ne_t, 'y':ne*1e-20, 'ylab':r'n$_e$ $10^{20}\,(m^{-3})$',\
               'lab':'proffit'}

print('Got ne from conf')

#te proffit
te_s = t.getNode(r'\tcv_shot::top.results.conf:te')
rho_conf_te = rho_s.data()
te = te_s.data()
te_t = te_s.getDimensionAt(1).data()
te0 = te[0,:]
te5 = te[20,:]
expdata[3] = {'nq':2, 'x':[te_t, te_t], \
              'y':[te0, te5], 'ylab':r'T$_e$ (keV)',\
              'lab':[r'$\rho_{POL}=0$', r'$\rho_{POL}=0.5$']}
expprof[3] = { 'x':rho_conf_te, 't':te_t, 'y':te*1e-3, 'ylab':r'T$_e$ (keV)$',\
               'lab':'proffit'}

print('Got te from conf')

# plot
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=10)
plt.rc('figure', facecolor='white')

ncol=3; nrow=3
times=np.linspace(0, np.max(expprof[0]['t']), ncol*nrow)
for kk in ['ne', 'te']:
    f, ax = plt.subplots(nrow, ncol, sharex='all', figsize=[4*ncol,3*nrow])
    if kk=='ne':
        f.suptitle('Ne for shot #'+str(shot), fontsize=20)
    else:
        f.suptitle('Te for shot #'+str(shot), fontsize=20)
    for i in range(ncol):
        for j in range(nrow):
            ind=i*nrow+j
            if kk=='ne':
                data1 = expprof[0]
                data2 = expprof[2]
                ax[i,j].set_ylabel(r'n$_e$ [10$^{20}$ m$^{-3}$]')
            else:
                data1 = expprof[1]
                data2 = expprof[3]
                ax[i,j].set_ylabel(r'T$_e$ [keV]')
                
            ttind1 = np.argmin(times[ind]-data1['t']>0.)
            ttind2 = np.argmin(times[ind]*1e-3-data2['t']>0.)
            
            ax[i,j].plot(data1['x'], data1['y'][:, ttind1], 'rx', lw=2.3, label='Thomson') #thomson
            ax[i,j].plot(data2['x'][ttind2,:], data2['y'][ttind2, :], 'b', lw=2.3, label='Conf') #conf
            ax[i,j].grid('on')
            ax[i,j].set_title('t={:2.2f} s'.format(times[ind]))
            
ax[0,0].set_xlim([0., 1.])
ax[0,0].legend(loc='lower left')
for j in range(ncol):
    ax[nrow-1,j].set_xlabel(r'$\rho_{POL}$');

plt.grid('on')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.show()
