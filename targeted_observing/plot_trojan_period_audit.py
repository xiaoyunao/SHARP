"""Scientific diagnostic figures for the November 11--13 catalog audit."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .audit_trojan_periods import fit
from .plot_style import apply_style, save_pair


def main():
    root=Path('targeted_observing/products/trojan_period_20261114')
    data=pd.read_csv(root/'differential_photometry.csv')
    data=data[(data.method=='PSF')&data.good]
    d=data[data.number==17365]
    t,y,e,n=[d[k].to_numpy() for k in ['hours','mag','err','night']]
    sensitivity=json.loads((root/'thymbraeus_sensitivity.json').read_text())
    period=sensitivity[0]['local_photometric_period_h']
    _,model,coef=fit(t,y,e,n,period,4)
    colours={20261111:'#1f77b4',20261112:'#ff7f0e',20261113:'#2ca02c'}
    markers={20261111:'o',20261112:'s',20261113:'^'}
    apply_style()
    # Visual styling only, based on the RAA 25745 panel in gotta_asteroid_1.
    # Preserve original panel content, phase range, fit order and sample filters.
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':1.3,
        'xtick.direction':'in','ytick.direction':'in','xtick.major.size':7,
        'ytick.major.size':7,'xtick.major.width':1.1,'ytick.major.width':1.1,
        'grid.linestyle':'--','grid.color':'.85','grid.linewidth':.9,
        'axes.labelsize':16,'axes.titlesize':17,'font.size':12})
    fig=plt.figure(figsize=(13.8,8.7),layout='constrained')
    gs=fig.add_gridspec(2,2,height_ratios=[.85,1.1])
    ax=fig.add_subplot(gs[0,:]);phase_ax=fig.add_subplot(gs[1,0]);scan=fig.add_subplot(gs[1,1])
    for night in np.unique(n):
        keep=n==night
        ax.errorbar(t[keep],y[keep]-np.median(y),e[keep],fmt=markers[night],
                    color=colours[night],ms=4,capsize=0,elinewidth=1.1,alpha=.85,label=f'Nov {str(night)[-2:]} (N={keep.sum()})')
        phase_ax.errorbar(t[keep]/period%1,y[keep]-np.median(y),e[keep],fmt=markers[night],
                         color=colours[night],ms=4,capsize=0,elinewidth=1.1,alpha=.85)
    ax.invert_yaxis();ax.set(xlabel='Hours since first observation (light-time corrected)',ylabel='Relative magnitude',
                             title='(17365) Thymbraeus: three observed nights')
    ax.legend(ncol=3,loc='lower left',fontsize=10,frameon=False);ax.grid(alpha=.55)
    ph=np.linspace(0,1,800);curve=np.full(len(ph),coef[0])
    for k in range(1,5):curve+=coef[2*k-1]*np.sin(2*np.pi*k*ph)+coef[2*k]*np.cos(2*np.pi*k*ph)
    phase_ax.plot(ph,curve-np.median(y),color='black',lw=2,zorder=0)
    phase_ax.invert_yaxis();phase_ax.set(xlabel='Brightness phase',ylabel='Relative magnitude',
        title=f'Brightness recurrence: {period:.3f} h')
    phase_ax.text(.5,.98,f'Double-peaked rotation candidate: {2*period:.2f} h',
                  transform=phase_ax.transAxes,va='top',ha='center',fontsize=9)
    phase_ax.grid(alpha=.55)
    grids=np.load(root/'period_grids.npz');p,c=grids['17365_PSF_H2_shared']
    keep=(p>=3)&(p<=40);scan.plot(p[keep],c[keep],color='#1f77b4',lw=1.3)
    for P,label in [(6.307,'6.3 h'),(12.641,'12.6 h'),(25.273,'25.3 h')]:
        cc=fit(t,y,e,n,P,2)[0];scan.plot(P,cc,'o',color='#ff8c00',ms=5)
        scan.annotate(label,(P,cc),xytext=(0,11),textcoords='offset points',ha='center',fontsize=10)
    scan.set(xlim=(3,40),ylim=(180,900),xlabel='Trial period (hours)',ylabel=r'$\chi^2$',
             title='Competing fits: two harmonics')
    scan.text(.97,.97,'25 h fit is unconstrained in phase gaps;\nits extrapolated amplitude is unphysical.',
              transform=scan.transAxes,ha='right',va='top',fontsize=9,
              bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=3))
    scan.grid(alpha=.55)
    save_pair(fig,root/'thymbraeus_period_diagnostic',
              'Exploratory differential catalog photometry. A 6.3 h brightness recurrence is not a unique spin period.')
    plt.close(fig)
    fig,axes=plt.subplots(5,1,figsize=(12,12),layout='constrained',sharex=True)
    names={17365:'Thymbraeus',2357:'Phereclos',51364:'2000 SU333',32513:'2001 OL31',4722:'Agelaos'}
    epoch=data.jd.min()
    for ax,number in zip(axes,[17365,2357,51364,32513,4722]):
        z=data[data.number==number]
        for night,g in z.groupby('night'):
            ax.errorbar((g.jd-epoch)*24,g.mag-z.mag.median(),g.err,
                        fmt=markers[night],color=colours[night],ms=3.5,elinewidth=1,alpha=.85)
        ax.invert_yaxis();ax.set_ylabel(r'$\Delta m$')
        ax.set_title(f'({number}) {names[number]}  |  {len(z)} accepted points',loc='left',fontsize=11,pad=5)
        ax.grid(alpha=.55)
    axes[-1].set_xlabel('Hours since first observation (UTC exposure midpoints)')
    save_pair(fig,root/'five_target_lightcurves','All five intended targets, PSF differential catalog photometry; gaps are real.')
    plt.close(fig)


if __name__=='__main__':main()
