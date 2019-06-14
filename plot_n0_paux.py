"""
routine to get paux from a shot
 and n0 from the output of get_n0LCS matlab routine
In Matlab:
n0 = get_n0LCS(<<SHOT>>, 1);
outdata=[n0.Time' n0.n0LCS];
save('/home/vallar/<<SHOT>>n0.dat', 'outdata', '-ascii')
"""

import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.interpolate as interp

shot=input('Shot number? ')
try:
    t=mds.Tree('tcv_shot', shot)
except:
    raise ValueError
expdata={}
col=['k', 'r','b']

#nbi power
pnbi_s = t.getNode(r'\ATLAS::NBH.DATA.MAIN_ADC:DATA')
pnbi = pnbi_s.data()[36,:]
pnbi_t = pnbi_s.getDimensionAt(0).data()
#ec power
pec_s = t.getNode(r'\results::toray.input:p_gyro')
pec = pec_s.data()[-1,:]*1e-3
pec_t = pec_s.getDimensionAt(0).data()
expdata[0] = {'nq':2, 'x':[pnbi_t, pec_t], \
              'y':[pnbi, pec], 'ylab':'',\
              'lab':['NBI (MW)', 'EC (MW)']}
print('Got paux')

try:
    d = np.loadtxt('/home/vallar/'+str(shot)+'n0.dat')
    print('Got n0')
except:
    print('No data for n0 found. Exit')
    raise ValueError
expdata[1] = {'nq':1, 'x':[d[:,0]], \
              'y':[d[:,1]*1e-16], 'ylab':r'$n_0^{LCS} [1/m^3]$',\
	      'lab':[r'$n_0^{LCS} [10^{16}/m^3]$']}

# plot
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('figure', facecolor='white')

f=plt.figure(figsize=(8,6))
ax = f.add_subplot(111)
f.suptitle('Shot #'+str(shot), fontsize=20)
data=expdata[0]
for el in range(data['nq']):
    ax.plot(data['x'][el], data['y'][el], col[el], lw=2.3, label=data['lab'][el])

data=expdata[1]
ax.plot(data['x'][0], data['y'][0], col[2], lw=2.3, label=data['lab'][0])
    
ax.set_xlabel(r'Time [s]'); ax.set_ylabel(r'')
ax.grid('on'); ax.legend(fontsize=20, loc='upper left')
ax.set_xlim([0., 2.4]); ax.set_ylim([0,6.])
plt.setp(ax.get_yticklabels()[0],visible=False)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.show()
