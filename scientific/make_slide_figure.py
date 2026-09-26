"""Large-label, single-row phase maps for the presentation."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from analyze import ROOT,FIG,COLORS,map_plot
from phaseguard import classify


def run():
    z=np.load(ROOT/'variational_map/variational_grid_N8.npz')
    uncertain=z['fidelity'][0]<.95;uncertain[0,:]=False
    hh,kk=np.where(uncertain)
    with plt.rc_context({'font.size':14,'axes.titlesize':20,'axes.labelsize':15,'xtick.labelsize':13,'ytick.labelsize':13}):
        fig,axs=plt.subplots(1,3,figsize=(13.5,5.2),layout='constrained')
        for ip,p in enumerate(z['p']):
            map_plot(axs[ip],z,classify(z['corr'][ip],float(z['threshold'])),f'p = {p:g}')
            axs[ip].scatter(z['kappa'][kk],z['h'][hh],s=30,facecolors='none',edgecolors='#303030',lw=1)
            axs[ip].set_xticks([0,.5,1]);axs[ip].set_yticks([0,.5,1,1.5,2])
        legend=[Line2D([0],[0],marker='s',linestyle='',color=color,markersize=11,label=label)
                for color,label in zip(COLORS,['Ferro-like','Antiphase-like','Weak order'])]
        legend.append(Line2D([0],[0],marker='o',linestyle='',color='#303030',markerfacecolor='none',
                             markersize=8,label='Clean fidelity < 0.95'))
        fig.legend(handles=legend,loc='outside lower center',ncols=4,fontsize=13)
        fig.savefig(FIG/'variational_maps_slides.png',dpi=180,bbox_inches='tight')
        plt.close(fig)


if __name__=='__main__':run()
