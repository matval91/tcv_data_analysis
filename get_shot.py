#!/usr/bin/env python2.7

# routine to get values from a shot and plot them
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.interpolate as interp
from matplotlib.gridspec import GridSpec
from utils.plot_utils import define_colors

colours, colours_old, styles, my_cmap, dpi = define_colors()

shot=input('Shot number? ')
time_req=input('Desired timeslice? ')

conn = mds.Connection('tcvdata.epfl.ch')
conn.openTree('tcv_shot', int(shot))

expdata={}
col=['k', 'r','b']
#plot 1
# ip; trapeze is measured, liuqe is calculated
sig=r'\magnetics::iplasma:trapeze'
iptr_s = conn.get(sig)
iptr = iptr_s.data()
iptr_t = conn.get('dim_of('+sig+',0)').data()
time = iptr_t
ind = np.argmin(time-time_req<0.)
sig=r'\results::i_p'
ipliu_s = conn.get(sig)
ipliu = ipliu_s.data()
ipliu_t = conn.get('dim_of('+sig+',0)').data()


expdata[0] = {'nq':2, 'x':[iptr_t, ipliu_t], \
              'y':[iptr*1e-3, ipliu*1e-3], 'ylab':r'$I_p$(kA)',\
              'lab':['trapeze', 'liuqe']}
print('Got Ip')

#plot 2
#vloop
sig=r'\magnetics::vloop'
vl_s = conn.get(sig)
vl = vl_s.data()[0,:]
vl_t = conn.get('dim_of('+sig+',0)').data()
vl_tf = np.linspace(min(vl_t), max(vl_t), 2000)
vl_p = interp.interp1d(vl_t, vl)
vl = vl_p(time)
vl_t = time
#btor
#toroidal vacuum field at geom. axis
sig=r'\magnetics::rbphi'
btor_s = conn.get(sig)
btor   = np.squeeze(btor_s.data()/0.88)
btor_t = conn.get('dim_of('+sig+',0)').data()
btor_par = interp.interp1d(btor_t, btor)

expdata[1] = {'nq':2, 'x':[vl_t, btor_t], \
              'y':[vl, btor], 'ylab':r'AU',\
              'lab':['V$_{loop}$ (V)', 'B$_{tor}(R_0)$ (T)'],
              'ylim':[-3, 5]}
print('Got vloop & Btor')

#plot 3
try:
    #nbi power
    sig=r'\ATLAS::NBH.DATA.MAIN_ADC:DATA'
    pnbi_s = conn.get(sig)
    pnbi = pnbi_s.data()[36,:]
    pnbi_t = conn.get('dim_of('+sig+',0)').data()
except:
    print('Most likely no beam is present')
    pnbi=np.zeros(10)
    pnbi_t=np.zeros(10)    
#ec power
try:
    sig=r'\results::toray.input:p_gyro'
    pec_s = conn.get(sig)
    pec = pec_s.data()[-1,:]*1e-3
    pec_t = conn.get('dim_of('+sig+',0)').data()
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
    sig=r'\base::pd:pd_001'
    da_s = conn.get(sig)
    da = da_s.data()
    da_t = conn.get('dim_of('+sig+',0)').data()
except:
    da=np.array([0])
    da_t=np.array([0])
#wmhd
sig=r'\results::total_energy:foo'
wmhd_s = conn.get(sig)
wmhd = wmhd_s.data()
wmhd_t = conn.get('dim_of('+sig+',0)').data()
expdata[3] = {'nq':2, 'x':[da_t, wmhd_t], \
              'y':[da*1000., wmhd], 'ylab':r'AU',\
              'lab':[r'D$_\alpha$', 'W$_{MHD}$']}
print('Got dalpha and wmhd')

#plot 5
#beta
sig=r'\results::beta_tor'
betat_s = conn.get(sig)
betat = betat_s.data()
betat_t = conn.get('dim_of('+sig+',0)').data()

#normalized beta = beta*a*btor/ip
# '100*beta/ip*1e6*b0*a_minor
sig=r'\results::r_max_psi'
rmax_s = conn.get(sig)
rmax = rmax_s.data()[:,-1]
sig=r'\results::r_min_psi'
rmin_s = conn.get(sig)
rmin = rmin_s.data()[:,-1]
print('Got stuff for betan')

a=(rmax-rmin)/2.
btor = np.abs(btor_par(betat_t))
a=a[~np.isnan(betat)]
btor=btor[~np.isnan(betat)]
ip_bt = np.abs(ipliu[~np.isnan(betat)])
betat_t = betat_t[~np.isnan(betat)]
betat = betat[~np.isnan(betat)]

betaN=100.*betat*a*1e6*btor/ip_bt
expdata[4] = {'nq':1, 'x':[betat_t], \
              'y':[betaN], 'ylab':r'$\beta_N$',\
              'lab':['']}
#plot 6
#ne
try:

    sig=r'\tcv_shot::top.results.conf:ne'
    ne_s=conn.get(sig)
    ne=ne_s.data()
    ne_t = conn.get('dim_of('+sig+',1)').data()
    ne_rho = conn.get('dim_of('+sig+',0)').data()
    ind_rho = np.argmin(ne_rho-0.5<0.)
    ne0 = ne[:,0]
    ne5 = ne[:,ind_rho]
    expdata[5] = {'nq':2, 'x':[ne_t, ne_t], \
                  'y':[ne0*1e-20, ne5*1e-20], 'ylab':r'n$_e$ $10^{20}\,(m^{-3})$',\
                  'lab':[r'$\rho=0$', r'$\rho=0.5$']}
    print('Got ne')
except:
    ne0=0;ne5=0.; ne_t=0.;
    expdata[5] = {'nq':2, 'x':[ne_t, ne_t], \
                  'y':[ne0*1e-20, ne5*1e-20], 'ylab':r'n$_e$ $10^{20}\,(m^{-3})$',\
                  'lab':[r'$\rho=0$', r'$\rho=0.5$']}
    print('No ne')

#plot 7
#te
try:
    sig=r'\tcv_shot::top.results.conf:te'
    te_s=conn.get(sig)
    te=te_s.data()
    te_t = conn.get('dim_of('+sig+',1)').data()
    te_rho = conn.get('dim_of('+sig+',0)').data()
    ind_rho = np.argmin(te_rho-0.5<0.)
    te0 = te[:,0]
    te5 = te[:,ind_rho]
    expdata[6] = {'nq':2, 'x':[te_t, te_t], \
                  'y':[te0, te5], 'ylab':r'T$_e$ (keV)',\
                  'lab':[r'$\rho=0$', r'$\rho=0.5$']}
    print('Got te')
except:
    te0=0;te5=0.; te_t=0.;
    expdata[6] = {'nq':2, 'x':[te_t, te_t], \
                  'y':[te0, te5], 'ylab':r'T$_e$ (keV)',\
                  'lab':[r'$\rho=0$', r'$\rho=0.5$']}
    print('No te')
#plot 8
#rz axis
sig=r'\results::z_axis'
zax_s = conn.get(sig)
zax = zax_s.data()*100.
zax_t = conn.get('dim_of('+sig+',0)').data()
sig=r'\results::r_axis'
rax_s = conn.get(sig)
rax = (rax_s.data()-0.88)*100.
rax_t = conn.get('dim_of('+sig+',0)').data()
expdata[7] = {'nq':2, 'x':[rax_t, zax_t], \
              'y':[rax, zax], 'ylab':r'Position (cm)',\
              'lab':[r'r$_{axis}$-0.88', r'z$_{axis}$']}
print('Got RZ axis')

#plot 9
#mhd
sig=r'\atlas::DT196_MHD_001:channel_067'
aaLFSz23_sect3=conn.get(sig)
sig=r'\atlas::DT196_MHD_001:channel_075'
aaLFSz23_sect11=conn.get(sig)
n1 = aaLFSz23_sect3
n1_sig = aaLFSz23_sect3.data() - aaLFSz23_sect11.data()
n2 = aaLFSz23_sect3
n2_sig = aaLFSz23_sect3.data() + aaLFSz23_sect11.data()
t_mhd=conn.get('dim_of('+sig+',0)').data()
mhd=dict()
mhd['n1']={};mhd['n2']={}
mhd['t']={}; mhd['tot']={}
mhd['n1']['sig'] = np.ravel(n1_sig)
mhd['n2']['sig'] = np.ravel(n2_sig)
mhd['t'] = np.ravel(t_mhd)
nfft=512
#[B,F,T]=specgram(gdat_data.data(:,i),nfft,1/mean(diff(gdat_data.t)),hanning(nfft),nfft/2);
for i in ['n1', 'n2']:
    mhd[i]['Pxx'], mhd[i]['freqs_1'], mhd[i]['bins'], mhd[i]['im'] = plt.specgram(mhd[i]['sig'], NFFT=nfft, Fs=1./np.mean(np.diff(t_mhd)))
    #plt.figure(); plt.specgram(mhd[i]['sig'], NFFT=nfft, Fs=1./np.mean(np.diff(t_mhd)))
mhd['tot']['Pxx'], mhd['tot']['freqs_1'], mhd['tot']['bins'], mhd['tot']['im'] = plt.specgram((mhd['n1']['sig']+mhd['n2']['sig'])/2., NFFT=nfft, Fs=1./np.mean(np.diff(t_mhd)))
expdata[9] = {}
expdata[9]['xbins'] = mhd['tot']['bins']
expdata[9]['ybins'] = mhd['tot']['freqs_1']
expdata[9]['values'] = mhd['tot']['Pxx']

#plt.figure(); plt.specgram((mhd['n1']['sig']+mhd['n2']['sig'])/2., NFFT=nfft, Fs=1./np.mean(np.diff(t_mhd)))
plt.close('all')
print('Got MHD')
#plot 10
#qpsi
sig=r'\results::q_psi'
qpsi_s = conn.get(sig)
qpsi = qpsi_s.data()
qpsi_t = conn.get('dim_of('+sig+',1)').data()
qpsi0 = qpsi[:,0]
qpsimin = np.min(qpsi, axis=1)
sig=r'\results::q_95'
qpsi95 = conn.get(sig).data()
expdata[8] = {'nq':3, 'x':[qpsi_t, qpsi_t, qpsi_t], \
              'y':[qpsi0, qpsimin, qpsi95], 'ylab':r'q',\
              'lab':[r'q$_0$', r'q$_{min}$', r'q$_{95}$'],\
              'ylim':[0, 10]}
print('Got q')

# RT signal to beam (need to connect to SCD tree and then RTC nodes)
# ask someone about it, it's not working.
# rt_s = t.getNode(r'\top::crpprt03:ethercat1:dac:dac_001')
# rt = rt_s.data()*10./32768.
# rt_t = rt_s conn.get('dim_of('+sig+',0)').data()
sig=r'\results::r_contour'
_R_s = conn.get(sig)
sig=r'\results::z_contour'
_Z_s = conn.get(sig)
_R = _R_s.data()
_z = _Z_s.data()
timeeq = conn.get('dim_of('+sig+',1)').data()
ind = np.argmin(timeeq-time_req<0)
R = _R[ind,:]; z=_z[ind,:]

_R = rax*0.01+0.88
_z = zax*0.01
timeeq = conn.get('dim_of('+sig+',1)')
ind = np.argmin(timeeq-time_req<0)
z0 = _z[ind]; R0=_R[ind]
rw=mds.Data.execute(r"static('r_v:in')").data()
zw=mds.Data.execute(r"static('z_v:in')").data()

if True:
    # plot
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=10)
    plt.rc('figure', facecolor='white')

    ncol=4; nrow=3
    #f, ax = plt.subplots(nrow, ncol, sharex='all', figsize=[4*ncol,3*nrow])
    f=plt.figure(figsize=[5*nrow,2*ncol])
    ax=np.empty((nrow, ncol), dtype='object')
    gs=GridSpec(nrow, ncol) 
    for i in range(nrow):
        for j in range(ncol-1):
            ax[i,j] = f.add_subplot(gs[i,j], sharex=ax[0,0])
    ax[0,ncol-1] = f.add_subplot(gs[0:2,ncol-1])
    ax[1,ncol-1] = []
    ax[2,ncol-1] = f.add_subplot(gs[2,ncol-1], sharex=ax[0,0])    
    
    f.suptitle('Shot #'+str(shot), fontsize=20)
    for i in range(nrow):
        for j in range(ncol):
            ind=i*nrow+j
            if len(expdata.keys())<=ind:
                break

            data=expdata[ind]
            if j==ncol-1:
                if i==0:
                    ax[i,j].plot(R,z, 'k');
                    ax[i,j].scatter(R0, z0, color='k' )
                    ax[i,j].plot(rw, zw, 'k', lw=2.3)
                    ax[i,j].plot(np.linspace(min(rw), max(rw), len(rw)), np.zeros(len(rw)), 'k--', lw=2.3);
                    ax[i,j].set_xlabel(r'R [m]'); ax[i,j].set_ylabel(r'z [m]')
                    ax[i,j].axis('equal')
                    continue
                elif i==1:
                    continue
                elif i==2:
                    ax[i,j].pcolormesh(t_mhd[0]+data['xbins'], data['ybins']*1e-3, np.log10(data['values']), cmap=my_cmap)
                    ax[i,j].set_ylabel(r'f [kHz]')
                    continue
            
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

    plt.grid('on')

    ax[0,0].set_xlim([0., 2.5])
    for j in range(ncol):
        ax[nrow-1,j].set_xlabel(r't (s)');

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

plt.show()
