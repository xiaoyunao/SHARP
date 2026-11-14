"""Redraw only the two requested overhead figures in the existing SHARP paper style."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullFormatter, StrMethodFormatter, AutoMinorLocator
from matplotlib import font_manager
import numpy as np

from .catalog import BASE

ROOT = BASE/'products/overhead_three_nights_20261114'
OUTPUT = ROOT/'paper_figures'
STYLE_PATH = BASE.parent/'paper_analysis_20260803/scripts/paper_style.py'


def load_style():
    spec = importlib.util.spec_from_file_location('sharp_existing_paper_style', STYLE_PATH)
    style = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(style)
    style.setup_style()
    plt.rcParams.update({'mathtext.fontset':'stix', 'axes.spines.top':True,
                         'axes.spines.right':True, 'axes.facecolor':'white',
                         'figure.facecolor':'white'})
    return style


def finish(ax, size=22):
    ax.set_axisbelow(True)
    ax.grid(which='major',alpha=.16,lw=.6)
    ax.tick_params(which='major',labelsize=size,direction='in',top=True,right=True)
    ax.tick_params(which='minor',direction='in',top=True,right=True,length=3.5,width=.8)


def draw():
    style=load_style();colors=[style.BLUE,style.ORANGE,style.GREEN]
    font=font_manager.findfont('Times New Roman',fallback_to_default=False)
    source=ROOT/'all_pairs.csv';rows=list(csv.DictReader(source.open()))
    nights=['20261111','20261112','20261113'];labels=['2026 Nov 11','2026 Nov 12','2026 Nov 13']
    summary=json.loads((ROOT/'summary.json').read_text())
    slews=[r for r in rows if r['kind']=='slew'];assert len(slews)==290
    OUTPUT.mkdir(exist_ok=True)

    fig,ax=plt.subplots(figsize=(10,8.6))
    ax.set(xscale='log',yscale='log',xlim=(9,310),ylim=(9,310))
    ax.plot([9,310],[9,310],color=style.INK,ls='--',lw=1.6,zorder=1)
    for night,color in zip(nights,colors):
        for context in ('within_block','block_boundary'):
            pp=[r for r in slews if r['night']==night and r['context']==context]
            x=np.array([float(r['theory_slew_s']) for r in pp]);y=np.array([float(r['inferred_slew_s']) for r in pp])
            assert np.all(x>0) and np.all(y>0)
            ax.scatter(x,y,s=32 if context=='within_block' else 85,
                marker='o' if context=='within_block' else '^',
                facecolors=color if context=='within_block' else 'none',
                edgecolors='none' if context=='within_block' else color,
                linewidths=1.5,alpha=.65 if context=='within_block' else 1,zorder=3)
    ticks=[10,20,30,50,100,200,300]
    for axis in (ax.xaxis,ax.yaxis):
        axis.set_major_locator(FixedLocator(ticks));axis.set_major_formatter(StrMethodFormatter('{x:g}'));axis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'Theoretical slew time [s]',labelpad=12)
    ax.set_ylabel(r'Measured slew time [s]',labelpad=12)
    handles=[Line2D([],[],ls='none',marker='o',color=c,ms=7,label=l) for c,l in zip(colors,labels)]
    handles += [Line2D([],[],ls='none',marker='^',mfc='none',mec=style.INK,ms=8,label='Block transition')]
    ax.legend(handles=handles,loc='upper left',fontsize=18,handletextpad=.45,borderaxespad=.75,labelspacing=.35)
    ax.text(175,205,'1:1',rotation=43,fontsize=19,color=style.INK)
    finish(ax)
    # The log main panel retains all points; a linear inset resolves the dense 30s region.
    inset=ax.inset_axes([.54,.10,.42,.37])
    inset.set(xlim=(24.5,34.5),ylim=(19.5,34.5))
    inset.plot([24.5,34.5],[24.5,34.5],'--',lw=1.1,color=style.INK)
    for night,color in zip(nights,colors):
        pp=[r for r in slews if r['night']==night and r['context']=='within_block']
        inset.scatter([float(r['theory_slew_s']) for r in pp],[float(r['inferred_slew_s']) for r in pp],s=17,color=color,alpha=.6,linewidths=0)
    inset.set_title('Neighbouring fields',fontsize=18,pad=8)
    inset.set_xticks([25,30]);inset.set_yticks([20,25,30]);finish(inset,16)
    for spine in inset.spines.values():spine.set_linewidth(1)
    fig.subplots_adjust(left=.16,right=.97,bottom=.17,top=.97)
    fig.text(.16,.035,r'Measured = inter-exposure gap $-$ 18 s; same-filter transfers only.',fontsize=14,color=style.INK)
    style.save_figure(fig,OUTPUT/'slew_theory_measured')

    # One contiguous chronological index: no artificial day-to-day gaps or jitter.
    repeats=sorted((r for r in rows if r['kind']=='repeat'),key=lambda r:(r['night'],int(r['to_sequence'])))
    assert len(repeats)==700
    fig,ax=plt.subplots(figsize=(13.2,6.5))
    ax.set(xlim=(.5,700.5),ylim=(16.4,23.2))
    start=1;medians={};boundaries=[];overflow=[]
    for night,color,label in zip(nights,colors,labels):
        rr=[r for r in repeats if r['night']==night];y=np.array([float(r['gap_s']) for r in rr])
        x=np.arange(start,start+len(rr));mid=float(np.median(y));medians[night]=mid
        assert abs(mid-summary['readout_per_night'][night]['median'])<1e-9
        good=y<=23.2
        ax.scatter(x[good],y[good],s=20,color=color,alpha=.7,linewidths=0,zorder=3)
        ax.hlines(mid,start-.5,start+len(rr)-.5,color=color,ls='--',lw=2.2,zorder=4)
        for ix,yy in zip(x[~good],y[~good]):
            ax.plot(ix,22.95,marker='^',ms=6,color=color,clip_on=False)
            overflow.append(dict(index=int(ix),gap_s=float(yy),night=night))
        middle=(x[0]+x[-1])/2
        if np.any(~good):
            label_high=', '.join(f'{v:.2f}' for v in y[~good])+' s'
            ax.text(middle,22.65,label_high,ha='center',va='center',fontsize=15,color=color)
        ax.text(middle,1.035,f'{label}\nMedian = {mid:.3f} s',
                transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=19,color=color)
        if start>1:
            ax.axvline(start-.5,color='0.6',lw=.8,ls=':',zorder=1);boundaries.append(start-.5)
        start+=len(rr)
    assert start==701 and boundaries==[266.5,462.5] and len(overflow)==3
    ax.set_xlabel('Consecutive measurement index',labelpad=10)
    ax.set_ylabel('CCD readout + control [s]',labelpad=12)
    ax.set_xticks(np.arange(0,701,100));ax.set_yticks([17,18,19,20,21,22,23]);ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    finish(ax)
    fig.subplots_adjust(left=.11,right=.98,bottom=.21,top=.83)
    fig.text(.11,.035,'The nights are concatenated without gaps. Off-scale intervals are marked at the top; all points enter the medians.',fontsize=13,color=style.INK)
    style.save_figure(fig,OUTPUT/'ccd_readout_sequence')
    verification=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        style_source=str(STYLE_PATH),style_source_sha256=hashlib.sha256(STYLE_PATH.read_bytes()).hexdigest(),font=font,
        slew_points=290,slew_axes='log-log; linear inset for neighbouring fields',
        readout_points=700,index_range=[1,700],night_boundaries=boundaries,medians_s=medians,
        off_scale_readout_points=overflow,analysis_unchanged=True)
    (OUTPUT/'plot_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification,indent=2))


if __name__=='__main__':draw()
