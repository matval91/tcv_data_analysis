import numpy as np
import MDSplus as mds
import utils.plot_utils as au
import matplotlib
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.interpolate as interp
col, col2, styles, my_cmap, dpi = au.define_colors()


#for now just need tor flux, pol flux, q, j, press

def get_tcv_eq(shot,time):
    """
    
    """
    conn = mds.Connection('tcvdata.epfl.ch')
    conn.openTree('tcv_shot', shot)
    
    #t = mds.Tree('tcv_shot', shot)
    #iptr_s = t.getNode(r'\magnetics::iplasma:trapeze')
    liuqe_flavour='LIUQE.M'
    _d={'label':'', 'x':[], 'y':[]}
    signals={'p':dict.copy(_d),
             'q':dict.copy(_d),
             'torflux':dict.copy(_d),
             'i':dict.copy(_d),
             'polflux':dict.copy(_d),
             'f':dict.copy(_d),
             'jbs':dict.copy(_d)}
    signals['p']['label']='p_rho'
    signals['q']['label']='q'
    signals['torflux']['label'] = 'tor_flux_tot'
    signals['polflux']['label'] = 'rho'
    signals['i']['label'] = 'i_p'
    signals['f']['label'] = 'rbtor_vac'
    signals['jbs']['label'] = r'\tcv_shot::top.results.ibs:jbsav'
    s_t=conn.get('dim_of(tcv_eq("i_p", "'+liuqe_flavour+'"))').data()
    ind = np.argmin(s_t-time<0.)    
    for i in signals.keys():
        if i=='i' or i=='f':
            _s = conn.get(r"tcv_eq('"+signals[i]['label']+"', '"+liuqe_flavour+"')")
            s = np.abs(_s.data()[ind])
            s_rho=s_t[ind]
        elif i=='polflux':
            _s = conn.get(r"tcv_eq('"+signals[i]['label']+"', '"+liuqe_flavour+"')")
            s = _s.data()
            s_rho = s
            _psiax = conn.get(r"tcv_eq('psi_axis', '"+liuqe_flavour+"')")
            psiax=_psiax.data()[ind]
            s = s**2.
            s = np.abs(s*psiax)/(2*np.pi)
        elif i=='jbs':
            t=mds.Tree('tcv_shot', shot)
            _s=t.getNode(signals[i]['label'])
            _ind = np.argmin(_s.getDimensionAt(1).data()-time<0.)    
            s = _s.data()[_ind,:]*-1.
            s_rho = _s.getDimensionAt(0).data()
            ps=interp.interp1d(s_rho,s)
            signals[i]['x'] = s_rho
            signals[i]['y'] = ps
            continue
        else:
            _s = conn.get(r"tcv_eq('"+signals[i]['label']+"', '"+liuqe_flavour+"')")
            s = _s.data()[ind,:]
            s_rho=conn.get('dim_of(tcv_eq("'+signals[i]['label']+'", "'+liuqe_flavour+'"),0)').data()

        signals[i]['x'] = s_rho
        signals[i]['y'] = s
    signals['torflux']['y']*=-1.
    signals['p']['y']*=1e-3
    return signals

def psi_to_phi(s):
    """
    Converts from poloidal to toroidal flux
    """
    locq=s['q']['y']
    locpsi=s['q']['x']**2
    
    if np.isinf(locq[-1]):
        locq[-1]=8
    phi = integrate.cumtrapz(locq,locpsi)
    phi = np.concatenate([[0], phi])
    phi=phi/np.max(phi)
    rhophi = phi**0.5
    for i in s:
        if i=='i' or i=='f':
            continue
        elif i=='jbs':
            s[i]['y'] = s[i]['y'](rhophi)
        #if i == 'q':
        #    s[i]['x'] == rhophi[:-1]; continue;
        s[i]['x']=rhophi
    return rhophi

def plot_tcv_eq(s, f=0):
    """
    """
    au.common_style()
    flag_leg=1
    ls='-'
    if not isinstance(f, matplotlib.figure.Figure):
        f = plt.figure(figsize=(12,10))
        axj = f.add_subplot(221)
        axp = f.add_subplot(222)
        axf = f.add_subplot(223)
        axc = f.add_subplot(224)
        axq = axf.twinx()  # instantiate a second axes that shares the same x-axis
        axff = axc.twinx()
    else:
        axj, axp, axf, axc, axq, axff =f.axes
        flag_leg=0
        ls='--'

    axj.plot(s['jbs']['x'], s['jbs']['y']*1e-3, 'r', label=r'BS', lw=2., linestyle=ls)
    
        
    axp.plot(s['p']['x'], s['p']['y'], 'k', label=r'Total exp', lw=2., linestyle=ls)
    axf.plot(s['torflux']['x'], s['torflux']['y'], 'b', label=r'Tor. Flux', lw=2., linestyle=ls)
    axf.plot(s['polflux']['x'], s['polflux']['y']*10., 'k', label=r'Pol. Flux x 10', lw=2., linestyle=ls)
    
    axq.plot(s['q']['x'], s['q']['y'], 'r', label=r'q', lw=2., linestyle=ls)
    axq.plot([0,1], [1,1], 'r--')
    axq.set_ylim([0,10]); axf.set_ylim([0,0.5])
    
    axc.plot(s['i']['x'], s['i']['y']*1e-3, 'ko', label=r'EXP CUR', lw=2., linestyle=ls)
    axff.plot(s['f']['x'], s['f']['y'], 'r^', label=r'EXP F', lw=2., linestyle=ls)
    axff.legend(loc='lower left')
    #========================================================
    # SET TICK LOCATION
    #========================================================
    if flag_leg==1:
        axq.tick_params(axis='y', labelcolor='r')
        axq.set_ylabel(r'q', color='r')
        axff.tick_params(axis='y', labelcolor='r')        
        au.limit_labels(axj, r'$\rho_{TOR}$', r'j [kA/m$^2$]','' )
        au.limit_labels(axp, r'$\rho_{TOR}$', r'p [kPa]','' )
        axf.set_xlabel(r'$\rho_{TOR}$')
        axf.set_ylabel(r'Fluxes (Wb/rad)')
        axj.legend(loc='best')
        axf.legend(loc='upper left'); axf.grid('on')
        axp.legend(loc='upper right'); axp.grid('on')
        axc.legend(loc='upper right'); axc.grid('on')
        au.limit_labels(axc, r'$t [s]$', r'I [kA]','' )        
        f.tight_layout()
        f.subplots_adjust(top=0.9)
    else:
        axq.set_yticklabels([])
 
        #if titles.get_text() == '':
            #f.suptitle('{} t={:.2f} s'.format(to.fname[-12:-4], time))
        #else:
            #newtitle = title+' {}'.format(to.fname[-12:-4])
        #    titles.set_text(newtitle)
    #except:
        #f.suptitle('{} t={:.2f} s'.format(to.fname[-12:-4], time))

    plt.show()  
