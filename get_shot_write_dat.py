# routine to get values from a shot and plot them
import MDSplus as mds
import numpy as np
import scipy.interpolate as interp

shot=input('Shot number? ')
t=mds.Tree('tcv_shot', shot)
expdata={}
dirt = '/home/vallar/'
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
print('Got Ip')

data = np.array([ipliu_t, ipliu])
np.savetxt(dirt+str(shot)+'ip.dat', data.T)
#plot 2
#vloop
vl_s = t.getNode(r'\magnetics::vloop')
vl = vl_s.data()[0,:]
vl_t = vl_s.getDimensionAt(0).data()
vl_tf = np.linspace(min(vl_t), max(vl_t), 2000)
vl_p = interp.interp1d(vl_t, vl)
vl = vl_p(time)
vl_t = time
data = np.array([vl_t, vl])
np.savetxt(dirt+str(shot)+'vl.dat', data.T)

#wmhd
wmhd_s = t.getNode(r'\results::total_energy:foo')
wmhd = wmhd_s.data()
wmhd_t = wmhd_s.getDimensionAt(0).data()
data = np.array([wmhd_t, wmhd])
np.savetxt(dirt+str(shot)+'wmhd.dat', data.T)                

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

btor_s = t.getNode(r'\magnetics::rbphi')
btor   = np.squeeze(btor_s.data()/0.88)
btor_t = btor_s.getDimensionAt(0).data()
btor_par = interp.interp1d(btor_t, btor)
a=(rmax-rmin)/2.
btor = btor_par(betat_t)
a=a[~np.isnan(betat)]
btor=btor[~np.isnan(betat)]
ip_bt = ipliu[~np.isnan(betat)]
betat_t = betat_t[~np.isnan(betat)]
betat = betat[~np.isnan(betat)]

betaN=100.*betat*a*1e6*btor/ip_bt

data=np.array([betat_t, betaN])
np.savetxt(dirt+str(shot)+'betaN.dat', data.T)


