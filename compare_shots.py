#!/usr/bin/env python2.7

# routine to get values from a shot and plot them
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.interpolate as interp

shot1 = 58832
shot2 = 58823

shot1=input('Shot number 1? ')
shot2=input('Shot number 2? ')

t1=mds.Tree('tcv_shot', shot1)
t2=mds.Tree('tcv_shot', shot2)

expdata1={}
expdata2={}
col=['k', 'r','b']
#plot 1
# ip; trapeze is measured, liuqe is calculated
iptr_s = t1.getNode(r'\magnetics::iplasma:trapeze')
iptr = iptr_s.data()
iptr_t = iptr_s.getDimensionAt(0).data()
time = iptr_t
ind = np.argmin(time-0.5<0.)
ipliu_s = t1.getNode(r'\results::i_p')
ipliu1 = ipliu_s.data()
ipliu_t = ipliu_s.getDimensionAt(0).data()

expdata1[0] = {'nq':1, 'x':[iptr_t], \
              'y':[-1.*iptr*1e-3], 'ylab':r'$I_p$(kA)',\
              'lab':['']}

# ip; trapeze is measured, liuqe is calculated
iptr_s = t2.getNode(r'\magnetics::iplasma:trapeze')
iptr = iptr_s.data()
iptr_t = iptr_s.getDimensionAt(0).data()
time = iptr_t
ind = np.argmin(time-0.5<0.)
ipliu_s = t2.getNode(r'\results::i_p')
ipliu2 = ipliu_s.data()
ipliu_t = ipliu_s.getDimensionAt(0).data()

expdata2[0] = {'nq':1, 'x':[iptr_t], \
              'y':[-1.*iptr*1e-3], 'ylab':r'$I_p$(kA)',\
              'lab':['']}
print('Got Ip')

#plot 2
#vloop
vl_s = t1.getNode(r'\magnetics::vloop')
vl = vl_s.data()[0,:]
vl_t = vl_s.getDimensionAt(0).data()
vl_tf = np.linspace(min(vl_t), max(vl_t), 200)
vl_p = interp.interp1d(vl_t, vl)
vl = vl_p(vl_tf)
vl_t = vl_tf

expdata1[1] = {'nq':1, 'x':[vl_t], \
              'y':[vl], 'ylab':r'',\
              'ylab':'V$_{loop}$ (V)',
              'ylim':[-3, 5]}
#vloop
vl_s = t2.getNode(r'\magnetics::vloop')
vl = vl_s.data()[0,:]
vl_t = vl_s.getDimensionAt(0).data()
vl_tf = np.linspace(min(vl_t), max(vl_t), 200)
vl_p = interp.interp1d(vl_t, vl)
vl = vl_p(vl_tf)
vl_t = vl_tf

expdata2[1] = {'nq':1, 'x':[vl_t], \
              'y':[vl], 'ylab':r'',\
              'ylab':'V$_{loop}$ (V)',
              'ylim':[-3, 5]}

print('Got vloop & Btor')

#plot 5
#beta
betat_s = t1.getNode(r'\results::beta_tor')
betat = betat_s.data()
betat_t = betat_s.getDimensionAt(0).data()

#btor
#toroidal vacuum field at geom. axis
btor_s = t1.getNode(r'\magnetics::rbphi')
btor   = np.squeeze(btor_s.data()/0.88)
btor_t = btor_s.getDimensionAt(0).data()
btor_par = interp.interp1d(btor_t, btor)

#normalized beta = beta*a*btor/ip
# '100*beta/ip*1e6*b0*a_minor
rmax_s = t1.getNode(r'\results::r_max_psi')
rmax = rmax_s.data()[:,-1]
rmin_s = t1.getNode(r'\results::r_min_psi')
rmin = rmin_s.data()[:,-1]
print('Got stuff for betan')

a=(rmax-rmin)/2.
btor = btor_par(betat_t)
a=a[~np.isnan(betat)]
btor=btor[~np.isnan(betat)]
ip_bt = ipliu1[~np.isnan(betat)]
betat_t = betat_t[~np.isnan(betat)]
betat = betat[~np.isnan(betat)]

betaN=100.*betat*a*1e6*btor/ip_bt
expdata1[2] = {'nq':1, 'x':[betat_t], \
              'y':[betaN], 'ylab':r'$\beta_N$',\
              'lab':['']}

#beta
betat_s = t2.getNode(r'\results::beta_tor')
betat = betat_s.data()
betat_t = betat_s.getDimensionAt(0).data()

#btor
#toroidal vacuum field at geom. axis
btor_s = t2.getNode(r'\magnetics::rbphi')
btor   = np.squeeze(btor_s.data()/0.88)
btor_t = btor_s.getDimensionAt(0).data()
btor_par = interp.interp1d(btor_t, btor)
#normalized beta = beta*a*btor/ip
# '100*beta/ip*1e6*b0*a_minor
rmax_s = t2.getNode(r'\results::r_max_psi')
rmax = rmax_s.data()[:,-1]
rmin_s = t2.getNode(r'\results::r_min_psi')
rmin = rmin_s.data()[:,-1]
print('Got stuff for betan')

a=(rmax-rmin)/2.
btor = btor_par(betat_t)
a=a[~np.isnan(betat)]
btor=btor[~np.isnan(betat)]
ip_bt = ipliu2[~np.isnan(betat)]
betat_t = betat_t[~np.isnan(betat)]
betat = betat[~np.isnan(betat)]

betaN=100.*betat*a*1e6*btor/ip_bt
expdata2[2] = {'nq':1, 'x':[betat_t], \
              'y':[betaN], 'ylab':r'$\beta_N$',\
              'lab':['']}

#plot 6
#ne
ne_s = t1.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:ne')
ne = ne_s.data()
ne_t = ne_s.getDimensionAt(0).data()
ne0 = ne[0,:]
ne5 = ne[20,:]
expdata1[3] = {'nq':1, 'x':[ne_t], \
              'y':[ne0*1e-20], 'ylab':r'n$_e$(0) $10^{20}\,(m^{-3})$',\
              'lab':[r'$\rho=0$']}
ne_s = t2.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:ne')
ne = ne_s.data()
ne_t = ne_s.getDimensionAt(0).data()
ne0 = ne[0,:]
ne5 = ne[20,:]
expdata2[3] = {'nq':1, 'x':[ne_t], \
              'y':[ne0*1e-20], 'ylab':r'n$_e$(0) $10^{20}\,(m^{-3})$',\
              'lab':[r'$\rho=0$']}

print('Got ne')

#plot 7
#te
te_s = t1.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:te')
te = te_s.data()
te_t = te_s.getDimensionAt(0).data()
te0 = te[0,:]*1e-3
te5 = te[20,:]*1e-3
expdata1[4] = {'nq':1, 'x':[te_t], \
              'y':[te0], 'ylab':r'T$_e$(0) (keV)',\
              'lab':[r'$\rho=0$']}
te_s = t2.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:te')
te = te_s.data()
te_t = te_s.getDimensionAt(0).data()
te0 = te[0,:]*1e-3
te5 = te[20,:]*1e-3
expdata2[4] = {'nq':1, 'x':[te_t], \
              'y':[te0], 'ylab':r'T$_e$(0) (keV)',\
              'lab':[r'$\rho=0$']}
print('Got te')

#plot 3
try:
    #nbi power
    pnbi_s = t1.getNode(r'\ATLAS::NBH.DATA.MAIN_ADC:DATA')
    pnbi = pnbi_s.data()[36,:]
    pnbi_t = pnbi_s.getDimensionAt(0).data()
except:
    pnbi=np.zeros(5)
    pnbi_t=np.zeros(5)    
try:
    #ec power
    pec_s = t1.getNode(r'\results::toray.input:p_gyro')
    pec = pec_s.data()[-1,:]*1e-3
    pec_t = pec_s.getDimensionAt(0).data()
except:
    pec=np.zeros((len(pnbi)))
    pec_t=np.zeros((len(pnbi)))
                
expdata1[5] = {'nq':2, 'x':[pnbi_t, pec_t], \
              'y':[pnbi, pec], 'ylab':r'P(MW)',\
              'lab':['NBI', 'EC']}

#plot 3

try:
    #nbi power
    pnbi_s = t2.getNode(r'\ATLAS::NBH.DATA.MAIN_ADC:DATA')
    pnbi = pnbi_s.data()[36,:]
    pnbi_t = pnbi_s.getDimensionAt(0).data()
except:
    pnbi=np.zeros(5)
    pnbi_t=np.zeros(5)    
#ec power
try:
    pec_s = t2.getNode(r'\results::toray.input:p_gyro')
    pec = pec_s.data()[-1,:]*1e-3
    pec_t = pec_s.getDimensionAt(0).data()
except:
    pec=np.zeros((len(pnbi)))
    pec_t=np.zeros((len(pnbi)))
                
expdata2[5] = {'nq':2, 'x':[pnbi_t, pec_t], \
              'y':[pnbi, pec], 'ylab':r'P(MW)',\
              'lab':['NBI', 'EC']}


# plot
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', facecolor='white')
stt=['-', '-']
ccc=['k','r']
ncol=3; nrow=2
f, ax = plt.subplots(nrow, ncol, sharex='all', figsize=[4*ncol,3*nrow])
f.text(0.3, 0.95,'Shot #'+str(shot1), fontsize=20, color='k')
f.text(0.6, 0.95,'Shot #'+str(shot2), fontsize=20, color='r')

for i in range(nrow):
    for j in range(ncol):
        ind=i+nrow*j
        if len(expdata1.keys())<=ind:
            break
        if ind==5:
            data = expdata1[ind];
            for el in range(data['nq']):
                col=['m','b']
                ax[i,j].plot(data['x'][el], data['y'][el], col[el], lw=1.5, ls=stile, label=data['lab'][el])
            ax[i,j].set_ylabel(r'P$_{AUX}$ [MW]')
            ax[i,j].grid('on')
            ax[i,j].legend(loc='best')
            if 'ylim' in data.keys():
                ax[i,j].set_ylim(data['ylim'])
            continue
        
        for dddi, ddd in enumerate([expdata1, expdata2]):
            data=ddd[ind]; stile = stt[dddi]
            col=['k','b']
            for el in range(data['nq']):
                ax[i,j].plot(data['x'][0], data['y'][0], ccc[dddi], lw=1.5)

        ax[i,j].set_ylabel(data['ylab'])
        ax[i,j].grid('on')
        if 'ylim' in data.keys():
            ax[i,j].set_ylim(data['ylim'])
            
        if data['nq']>1:
            ax[i,j].legend(loc='best')

            
ax[0,0].set_xlim([0., 3])
for j in range(ncol):
    ax[nrow-1,j].set_xlabel(r't (s)');

plt.grid('on')
plt.tight_layout()
plt.subplots_adjust(top=0.9)

_R_s = t1.getNode(r'\results::r_contour')
_Z_s = t1.getNode(r'\results::z_contour')
_R1 = _R_s.data()
_z1 = _Z_s.data()
t = _R_s.getDimensionAt(1).data()
ind = np.argmin(t-0.8<0)
R1 = _R1[ind,:]; z1=_z1[ind,:]
_R = t1.getNode(r'\results::r_axis')
_z = t1.getNode(r'\results::z_axis')
_z1 = _z.data(); _R1=_R.data()
t = _z.getDimensionAt(0)
ind = np.argmin(t-0.8<0)
z01 = _z1[ind]; R01=_R1[ind]

_R_s = t2.getNode(r'\results::r_contour')
_Z_s = t2.getNode(r'\results::z_contour')
_R2 = _R_s.data()
_z2 = _Z_s.data()
t = _R_s.getDimensionAt(1).data()
ind = np.argmin(t-1.25<0)
R2 = _R2[ind,:]; z2=_z2[ind,:]
_R = t2.getNode(r'\results::r_axis')
_z = t2.getNode(r'\results::z_axis')
_z2 = _z.data(); _R2=_R.data()
t = _z.getDimensionAt(0)
ind = np.argmin(t-1.25<0)
z02 = _z2[ind]; R02=_R2[ind]

#wall = np.loadtxt('/home/vallar/TCV_wall/TCV_FW_coord.dat')
rw=mds.Data.execute(r"static('r_v:in')").data()
zw=mds.Data.execute(r"static('z_v:in')").data()
f=plt.figure(figsize=(3,6)); ax=f.add_subplot(111)
ax.plot(R1,z1, 'r', label=str(shot1)); ax.plot(R2,z2,'k', label=str(shot2))
ax.scatter(R01, z01, color='r', ); ax.scatter(R02, z02, color='k')
ax.plot(rw, zw, 'k', lw=2.3)
ax.plot(np.linspace(min(rw), max(rw), len(rw)), np.zeros(len(rw)), 'k--', lw=2.3);
ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
ax.axis('equal')
ax.legend(loc='best')
f.tight_layout()
plt.show()
