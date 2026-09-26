"""Analyze the independent 15x15 short-circuit map after it is copied here."""
import json
import numpy as np
import matplotlib.pyplot as plt
from analyze import ROOT,FIG,COLORS,CMAP,map_plot,reference_labels,Line2D
from phaseguard import order,classify,crossing


def run():
    folder=ROOT/'variational_map'
    z=np.load(folder/'variational_grid_N8.npz')
    assert z['complete'].all() and z['corr'].shape==(3,15,15,8)
    t=float(z['threshold']);labels=np.array([classify(c,t) for c in z['corr']])
    exact=classify(z['exact_corr'],t)
    f,a,_=order(z['corr'])
    uncertain=z['fidelity'][0]<.95
    uncertain[0,:]=False  # h=0 needs a ground-space metric, not one state.
    uh,uk=np.where(uncertain)
    def mark_uncertain(ax):
        ax.scatter(z['kappa'][uk],z['h'][uh],s=20,facecolors='none',edgecolors='#303030',lw=.7)
    ref,mask=reference_labels(z['kappa'],z['h'])
    result=dict(grid=[15,15],cnot_count=28,threshold=t,
                clean_label_agreement_with_exact=float(np.mean(labels[0]==exact)),
                clean_reference_agreement_outside_floating=float(np.mean(labels[0][mask]==ref[mask])),
                mean_clean_energy_error_per_site=float(z['clean_energy_error_per_site'].mean()),
                max_clean_energy_error_per_site=float(z['clean_energy_error_per_site'].max()),
                mean_clean_fidelity_h_positive=float(z['fidelity'][0,1:].mean()),
                min_clean_fidelity_h_positive=float(z['fidelity'][0,1:].min()),
                optimizer_success_count=int(z['optimizer_success'].sum()),
                max_trace_error=float(abs(z['trace']-1).max()),
                high_energy_error_count=int(np.sum(z['clean_energy_error_per_site']>.0025)),
                low_clean_fidelity_positive_field_count=int(uncertain.sum()),
                noise=[],boundaries=[])
    for ip,p in enumerate(z['p']):
        result['noise'].append(dict(p=float(p),counts=np.bincount(labels[ip].ravel(),minlength=3).tolist(),
                                   clean_VQE_label_agreement=float(np.mean(labels[ip]==labels[0])),
                                   exact_clean_label_agreement=float(np.mean(labels[ip]==exact)),
                                   ferro_recall=float(np.mean(labels[ip][labels[0]==0]==0)),
                                   anti_recall=float(np.mean(labels[ip][labels[0]==1]==1)),
                                   mean_fidelity_h_positive=float(z['fidelity'][ip,1:].mean()),
                                   mean_purity=float(z['purity'][ip].mean())))
    for target in [0,.3,.6,.8,1.]:
        j=int(np.argmin(abs(z['kappa']-target)));k=float(z['kappa'][j])
        for ip,p in enumerate(z['p']):
            val=f[ip,:,j] if k<.5 else a[ip,:,j]
            value=crossing(z['h'],val,t)
            bracket=None if value is None else min(int(np.searchsorted(z['h'],value)),len(z['h'])-1)
            flag=None if bracket is None else bool(uncertain[max(0,bracket-1):bracket+1,j].any())
            result['boundaries'].append(dict(kappa=k,p=float(p),h=value,uncertain_preparation_in_bracket=flag))
    (ROOT/'data/variational_summary.json').write_text(json.dumps(result,indent=2))
    legend=[Line2D([0],[0],marker='s',linestyle='',color=col,markersize=9,label=name)
            for col,name in zip(COLORS,['Ferro-like','Antiphase-like','Weak order / para-like'])]
    legend.append(Line2D([0],[0],marker='o',linestyle='',color='#303030',markerfacecolor='none',
                         markersize=6,label='Clean fidelity < 0.95'))
    fig,axs=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    for ip,p in enumerate(z['p']):
        map_plot(axs[0,ip],z,labels[ip],f'p = {p:g}')
        mark_uncertain(axs[0,ip])
        im=map_plot(axs[1,ip],z,np.maximum(f[ip],a[ip]),'Order signal',False)
    fig.colorbar(im,ax=axs[1,:],label='max(F, A)',shrink=.8)
    fig.legend(handles=legend,loc='outside lower center',ncols=4)
    fig.suptitle('Short preparation keeps more order | N = 8 · 28 CNOTs',fontsize=17,fontweight='bold')
    fig.savefig(FIG/'variational_triptych.png',dpi=170,bbox_inches='tight');plt.close(fig)
    for ip,p in enumerate(z['p']):
        fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
        map_plot(axs[0],z,labels[ip],f'Fixed order classifier · p = {p:g}')
        mark_uncertain(axs[0])
        im=map_plot(axs[1],z,np.maximum(f[ip],a[ip]),'Order signal max(F, A)',False)
        fig.colorbar(im,ax=axs[1],label='Normalized order signal',shrink=.85)
        axs[0].legend(handles=legend,fontsize=8,loc='upper left')
        fig.suptitle('PhaseGuard | N = 8 · 28 noisy preparation CNOTs',fontsize=15,fontweight='bold')
        fig.text(.5,-.015,'15 × 15 grid · fixed clean VQE parameters · solid: Ising/KT; dashed: PT reference',ha='center',fontsize=9)
        fig.savefig(FIG/f'vqe_phase_p{p:.2f}.png',dpi=170,bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,target in zip(axs,[.3,.8]):
        j=int(np.argmin(abs(z['kappa']-target)));k=float(z['kappa'][j])
        for ip,p in enumerate(z['p']):
            y=f[ip,:,j] if k<.5 else a[ip,:,j]
            ax.plot(z['h'],y,'-o',ms=3,label=f'p = {p:g}')
        ax.axhline(t,color='k',ls='--',lw=.8,label='Fixed threshold')
        ax.set(xlabel='Transverse field h',ylabel='Order signal',title=f'28-CNOT preparation · κ={k:.3f}',xlim=(0,1.3))
        ax.legend(fontsize=8)
    fig.savefig(FIG/'variational_boundaries.png',dpi=180,bbox_inches='tight');plt.close(fig)
    # Keep optimizer error separate from simulated gate noise.
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    im=axs[0].pcolormesh(z['kappa'],z['h'],z['clean_energy_error_per_site'],shading='nearest',cmap='magma')
    fig.colorbar(im,ax=axs[0],label='(EVQE − Eexact)/N')
    im=axs[1].pcolormesh(z['kappa'],z['h'][1:],z['fidelity'][0,1:],shading='nearest',cmap='viridis',vmin=0,vmax=1)
    fig.colorbar(im,ax=axs[1],label='Clean fidelity; h=0 excluded')
    for ax in axs:ax.set(xlabel='Frustration κ',ylabel='Transverse field h')
    axs[0].set_title('Preparation error before noise');axs[1].set_title('State check before noise')
    fig.savefig(FIG/'variational_quality.png',dpi=170,bbox_inches='tight');plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
