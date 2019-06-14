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
# ip; trapeze is measured, liuqe is calculated
iptr_s = t.getNode(r'\magnetics::iplasma:trapeze')
iptr = iptr_s.data()
iptr_t = iptr_s.getDimensionAt(0).data()
time = iptr_t
ind = np.argmin(time-0.5<0.)
ipliu_s = t.getNode(r'\results::i_p')
ipliu = ipliu_s.data()
ipliu_t = ipliu_s.getDimensionAt(0).data()


expdata[0] = {'nq':2, 'x':[iptr_t, ipliu_t], \
              'y':[iptr*1e-3, ipliu*1e-3], 'ylab':r'$I_p$(kA)',\
              'lab':['trapeze', 'liuqe']}
print('Got Ip')

#plot 2
#vloop
vl_s = t.getNode(r'\magnetics::vloop')
vl = vl_s.data()[0,:]
vl_t = vl_s.getDimensionAt(0).data()
vl_tf = np.linspace(min(vl_t), max(vl_t), 2000)
vl_p = interp.interp1d(vl_t, vl)
vl = vl_p(time)
vl_t = time
#btor
#toroidal vacuum field at geom. axis
btor_s = t.getNode(r'\magnetics::rbphi')
btor   = np.squeeze(btor_s.data()/0.88)
btor_t = btor_s.getDimensionAt(0).data()
btor_par = interp.interp1d(btor_t, btor)

expdata[1] = {'nq':2, 'x':[vl_t, btor_t], \
              'y':[vl, btor], 'ylab':r'AU',\
              'lab':['V$_{loop}$ (V)', 'B$_{tor}(R_0)$ (T)'],
              'ylim':[-3, 5]}
print('Got vloop & Btor')

#plot 3
#nbi power
pnbi_s = t.getNode(r'\ATLAS::NBH.DATA.MAIN_ADC:DATA')
pnbi = pnbi_s.data()[36,:]
pnbi_t = pnbi_s.getDimensionAt(0).data()
#ec power
try:
    pec_s = t.getNode(r'\results::toray.input:p_gyro')
    pec = pec_s.data()[-1,:]*1e-3
    pec_t = pec_s.getDimensionAt(0).data()
except:
    pec=np.zeros((len(pnbi)))
    pec_t=np.zeros((len(pnbi)))
                
expdata[2] = {'nq':2, 'x':[pnbi_t, pec_t], \
              'y':[pnbi, pec], 'ylab':r'P(MW)',\
              'lab':['NBI', 'EC']}
print('Got paux')

#plot 4
#d-alpha
try:
    da_s = t.getNode(r'\base::pd:pd_001')
    da = da_s.data()
    da_t = da_s.getDimensionAt(0).data()
except:
    da=np.array([0])
    da_t=np.array([0])
#wmhd
wmhd_s = t.getNode(r'\results::total_energy:foo')
wmhd = wmhd_s.data()
wmhd_t = wmhd_s.getDimensionAt(0).data()
expdata[3] = {'nq':2, 'x':[da_t, wmhd_t], \
              'y':[da*1000., wmhd], 'ylab':r'AU',\
              'lab':[r'D$_\alpha$', 'W$_{MHD}$']}
print('Got dalpha and wmhd')

#plot 5
#beta
betat_s = t.getNode(r'\results::beta_tor')
betat = betat_s.data()
betat_t = betat_s.getDimensionAt(0).data()

#normalized beta = beta*a*btor/ip
# '100*beta/ip*1e6*b0*a_minor
rmax_s = t.getNode(r'\results::r_max_psi')
rmax = rmax_s.data()[:,-1]
rmin_s = t.getNode(r'\results::r_min_psi')
rmin = rmin_s.data()[:,-1]
print('Got stuff for betan')

a=(rmax-rmin)/2.
btor = btor_par(betat_t)
a=a[~np.isnan(betat)]
btor=btor[~np.isnan(betat)]
ip_bt = ipliu[~np.isnan(betat)]
betat_t = betat_t[~np.isnan(betat)]
betat = betat[~np.isnan(betat)]

betaN=100.*betat*a*1e6*btor/ip_bt
expdata[4] = {'nq':1, 'x':[betat_t], \
              'y':[betaN], 'ylab':r'$\beta_N$',\
              'lab':['']}
#plot 6
#ne
ne_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:ne')
ne = ne_s.data()
ne_t = ne_s.getDimensionAt(0).data()
ne0 = ne[0,:]
ne5 = ne[20,:]
expdata[5] = {'nq':2, 'x':[ne_t, ne_t], \
              'y':[ne0*1e-20, ne5*1e-20], 'ylab':r'n$_e$ $10^{20}\,(m^{-3})$',\
              'lab':[r'$\rho=0$', r'$\rho=0.5$']}
print('Got ne')

#plot 7
#te
te_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:te')
te = te_s.data()
te_t = te_s.getDimensionAt(0).data()
te0 = te[0,:]
te5 = te[20,:]
expdata[6] = {'nq':2, 'x':[te_t, te_t], \
              'y':[te0, te5], 'ylab':r'T$_e$ (keV)',\
              'lab':[r'$\rho=0$', r'$\rho=0.5$']}
print('Got te')

#plot 8
#rz axis
zax_s = t.getNode(r'\results::z_axis')
zax = zax_s.data()
zax_t = zax_s.getDimensionAt(0)
rax_s = t.getNode(r'\results::r_axis')
rax = rax_s.data()-0.88
rax_t = rax_s.getDimensionAt(0)
expdata[7] = {'nq':2, 'x':[rax_t, zax_t], \
              'y':[rax, zax], 'ylab':r'Position (m)',\
              'lab':[r'r$_{axis}$-0.88', r'z$_{axis}$']}
print('Got RZ axis')

#plot 9
#psi
psi=t.getNode(r'\results::psi')


#plot 10
#qpsi
qpsi_s = t.getNode(r'\results::q_psi')
qpsi = qpsi_s.data()
qpsi_t = qpsi_s.getDimensionAt(1).data()
qpsi0 = qpsi[:,0]
qpsimin = np.min(qpsi, axis=1)
qpsi95 = t.getNode(r'\results::q_95').data()
expdata[8] = {'nq':3, 'x':[qpsi_t, qpsi_t, qpsi_t], \
              'y':[qpsi0, qpsimin, qpsi95], 'ylab':r'q',\
              'lab':[r'q$_0$', r'q$_{min}$', r'q$_{95}$'],\
              'ylim':[0, 10]}
print('Got q')





# RT signal to beam (need to connect to SCD tree and then RTC nodes)
# ask someone about it, it's not working.
# rt_s = t.getNode(r'\top::crpprt03:ethercat1:dac:dac_001')
# rt = rt_s.data()*10./32768.
# rt_t = rt_s.getDimensionAt(0).data()



# plot
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=10)
plt.rc('figure', facecolor='white')

ncol=3; nrow=3
f, ax = plt.subplots(nrow, ncol, sharex='all', figsize=[4*ncol,3*nrow])
f.suptitle('Shot #'+str(shot), fontsize=20)
for i in range(ncol):
    for j in range(nrow):
        ind=i*nrow+j
        if len(expdata.keys())<=ind:
            break

        data=expdata[ind]
        for el in range(data['nq']):
            ax[i,j].plot(data['x'][el], data['y'][el], col[el], lw=1.5, label=data['lab'][el])
        if ind==7: #plot 0 position with raxis and zaxis
            ax[i,j].plot(data['x'][0], np.zeros(len(data['x'][0])), 'k--')
        if ind==8: #plot q=1
            ax[i,j].plot(data['x'][0], np.ones(len(data['x'][0])), 'k--')
        ax[i,j].set_ylabel(data['ylab'])
        ax[i,j].grid('on')
        if 'ylim' in data.keys():
            ax[i,j].set_ylim(data['ylim'])
            
        if data['nq']>1:
            ax[i,j].legend(loc='best')

            
ax[0,0].set_xlim([0., 2.5])
for j in range(ncol):
    ax[nrow-1,j].set_xlabel(r't (s)');

plt.grid('on')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.show()
