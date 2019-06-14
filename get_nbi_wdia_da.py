# routine to get pnbi, wmhd and dalpha from a shot
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.interpolate as interp
shot=input('Shot number? ')
t=mds.Tree('tcv_shot', shot)

#nbi power
pnbi_s = t.getNode(r'\ATLAS::NBH.DATA.MAIN_ADC:DATA')
pnbi = pnbi_s.data()[36,:]
pnbi_t = pnbi_s.getDimensionAt(0).data()

#ec power
pec_s = t.getNode(r'\results::toray.input:p_gyro')
pec = pec_s.data()[-1,:]*1e-3
pec_t = pec_s.getDimensionAt(0).data()

#d-alpha
da_s = t.getNode(r'\base::pd:pd_001')
da = da_s.data()
da_t = da_s.getDimensionAt(0).data()

#vloop
vl_s = t.getNode(r'\magnetics::vloop')
vl = vl_s.data()[0,:]
vl_t = vl_s.getDimensionAt(0).data()
vl_tf = np.linspace(min(vl_t), max(vl_t), 2000)
vl_p = interp.interp1d(vl_t, vl)
vl = vl_p(vl_tf)
vl_t = vl_tf
#wmhd
wmhd_s = t.getNode(r'\results::total_energy:foo')
wmhd = wmhd_s.data()
wmhd_t = wmhd_s.getDimensionAt(0).data()

#beta
betat_s = t.getNode(r'\results::beta_tor')
betat = betat_s.data()
betat_t = betat_s.getDimensionAt(0).data()

# RT signal to beam (need to connect to SCD tree and then RTC nodes)
# ask someone about it, it's not working.
# rt_s = t.getNode(r'\top::crpprt03:ethercat1:dac:dac_001')
# rt = rt_s.data()*10./32768.
# rt_t = rt_s.getDimensionAt(0).data()

#normalized beta = beta*a*btor/ip
# '100*beta/ip*1e6*b0*a_minor
rmax_s = t.getNode(r'\results::r_max_psi')
rmax = rmax_s.data()[:,-1]
rmin_s = t.getNode(r'\results::r_min_psi')
rmin = rmin_s.data()[:,-1]
a=(rmax-rmin)/2.

#toroidal vacuum field at geom. axis
btor_s = t.getNode(r'\magnetics::rbphi')
btor   = np.squeeze(btor_s.data()/0.88)
btor_t = btor_s.getDimensionAt(0).data()
btor_par = interp.interp1d(btor_t, btor)
btor = btor_par(betat_t)
#ip
ip_s = t.getNode(r'\results::i_p')
ip = ip_s.data()

a=a[~np.isnan(betat)]
btor=btor[~np.isnan(betat)]
ip = ip[~np.isnan(betat)]
betat_t = betat_t[~np.isnan(betat)]
betat = betat[~np.isnan(betat)]

betaN=100.*betat*a*1e6*btor/ip

# plot
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('figure', facecolor='white')

f=plt.figure(); ax=f.add_subplot(111)
#ax.set_title(r'Shot #'+str(shot)+' - NBH L-H transition, ELMy')
ax.set_title(r'Shot #'+str(shot)+' - RT beam control')
ax.plot(da_t, da*0.1, 'b', lw=2., label=r'D$_\alpha$ (AU)')
#ax.plot(vl_t, vl*0.2, 'b', lw=2., label=r'V$_{loop}$ (AU)')
#ax.plot(pec_t, pec, 'g--', lw=2., label=r'EC (MW)')
ax.plot(pnbi_t, pnbi, 'r', lw=2.3, label=r'NBI (MW)')
#ax.plot(wmhd_t, wmhd*1e-3, 'k', lw=2., label=r'W$_{DML}$ (kJ)')
ax.plot(betat_t, betaN, 'k', lw=2., label=r'$\beta_N$')

ax.set_xlim([0.6, 1.3]); ax.set_ylim([0, 2])
ax.set_xlabel(r't (s)'); ax.set_ylabel('AU')
plt.grid('on')

# Create your ticker object with M ticks
M = 4
yticks = ticker.MaxNLocator(M)
xticks = ticker.MaxNLocator(M)
# tick positions with set_minor_locator.
ax.yaxis.set_major_locator(yticks)
ax.xaxis.set_major_locator(xticks)
#==============================================
#Removing first point of y-axis
plt.setp(ax.get_yticklabels()[0], visible=False) 

plt.tight_layout()
plt.legend(loc='best', fontsize='large')
plt.show()
