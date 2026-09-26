"""Build every figure and numeric summary from the saved experiment data."""
import csv
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from phaseguard import ROOT,ANNNI,reference,attenuation,order,classify,crossing

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.facecolor":"#fcfcfb","axes.facecolor":"#fcfcfb",
                     "savefig.facecolor":"#fcfcfb","axes.titleweight":"bold"})
COLORS=["#267c83","#e99b44","#dce2e9"]
CMAP=ListedColormap(COLORS)
FIG=ROOT/"figures"; FIG.mkdir(exist_ok=True)


def threshold(n):
    m=ANNNI(n); _,s,_=m.ground(0,1)
    return float(order(m.observables(s)[0])[0])


def overlay(ax):
    k=np.linspace(0,1,401); a,b,c=reference(k)
    ax.plot(k,a,"k-",lw=1.1)
    ax.plot(k,c,"k-",lw=1.1)
    ax.plot(k,b,"k--",lw=1.1)
    ax.plot(.5,0,"ko",ms=4)
    ax.set(xlim=(0,1),ylim=(0,2),xlabel="Frustration κ",ylabel="Transverse field h")


def map_plot(ax,grid,values,title,categorical=True):
    if categorical:
        im=ax.pcolormesh(grid["kappa"],grid["h"],values,cmap=CMAP,vmin=-.5,vmax=2.5,shading="nearest")
    else:
        im=ax.pcolormesh(grid["kappa"],grid["h"],values,cmap="magma",vmin=0,vmax=1,shading="nearest")
    overlay(ax); ax.set_title(title)
    return im


def reference_labels(ks,hs):
    k,h=np.meshgrid(ks,hs); left,lower,upper=reference(ks)
    labels=np.full(k.shape,2)
    labels[(k<.5)&(h<left[None,:])]=0
    labels[(k>.5)&(h<upper[None,:])]=1
    mask=(k!=.5)&(h>0)&~((k>.5)&(h>=lower[None,:])&(h<=upper[None,:]))
    return labels,mask


def run():
    clean=np.load(ROOT/"data/grid_N12.npz")
    prep=np.load(ROOT/"data/preparation_grid_N8.npz")
    mit=np.load(ROOT/"data/mitigation_grid_N8.npz")
    t8,t12=threshold(8),threshold(12)
    c=prep["corr"]
    labels=np.array([classify(x,t8) for x in c])
    f,a,s=order(c)
    zne=np.array([2*c[1]-mit["corr"][1],2*c[2]-mit["corr"][2]])
    zlabels=np.array([classify(x,t8) for x in zne])
    ref,mask=reference_labels(prep["kappa"],prep["h"])
    counts=[{name:int(np.sum(l==i)) for i,name in enumerate(["ferro","anti","weak_order"])} for l in labels]
    summary={"N8_threshold":t8,"N12_threshold":t12,"N8_grid_shape":list(labels.shape[1:]),
             "N12_grid_shape":list(clean["energy"].shape),"preparation_CNOT_range":[int(prep["cnot_count"].min()),int(prep["cnot_count"].max())],
             "preparation_max_trace_error":float(abs(prep["trace"]-1).max()),
             "N12_max_residual":float(clean["residual"].max()),"phase_counts":counts,
             "reference_agreement_N8_clean_outside_floating_band":float(np.mean(labels[0][mask]==ref[mask])),
             "grid_spacing":{"N8_h":float(prep["h"][1]-prep["h"][0]),"N12_h":float(clean["h"][1]-clean["h"][0])},
             "noise":[],"threshold_sensitivity":[]}
    for ip,p in enumerate([.01,.05],1):
        record={"p":p,"label_agreement_with_clean":float(np.mean(labels[ip]==labels[0])),
                "ferro_recall":float(np.mean(labels[ip][labels[0]==0]==0)),
                "anti_recall":float(np.mean(labels[ip][labels[0]==1]==1)),
                "correlation_MAE":float(np.mean(abs(c[ip,...,1:]-c[0,...,1:]))),
                "zne_label_agreement_with_clean":float(np.mean(zlabels[ip-1]==labels[0])),
                "zne_recovered_ordered_labels":int(np.sum((zlabels[ip-1]==labels[0])&(labels[0]!=2))),
                "zne_false_positive_ordered_labels":int(np.sum((zlabels[ip-1]!=2)&(labels[0]==2))),
                "zne_correlation_MAE":float(np.mean(abs(zne[ip-1,...,1:]-c[0,...,1:]))),
                "mean_ground_state_fidelity":float(prep["fidelity"][ip].mean()),
                "mean_purity":float(prep["purity"][ip].mean())}
        summary["noise"].append(record)
    for t in [.25,.4,t8]:
        summary["threshold_sensitivity"].append({"threshold":t,"ordered_points":[int(np.sum(classify(x,t)!=2)) for x in c]})
    legend=[Line2D([0],[0],marker="s",linestyle="",color=col,markersize=10,label=name)
            for col,name in zip(COLORS,["Ferro-like","Antiphase-like","Weak order / para-like"])]
    for ip,p in enumerate(prep["p"]):
        fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
        map_plot(axs[0],prep,labels[ip],f"Fixed order classifier · p = {p:g}")
        im=map_plot(axs[1],prep,np.maximum(f[ip],a[ip]),"Order signal max(F, A)",False)
        fig.colorbar(im,ax=axs[1],label="Normalized order signal",shrink=.85)
        axs[0].legend(handles=legend,fontsize=8,loc="upper left")
        fig.suptitle("PhaseGuard | N = 8 · noise after every preparation CNOT",fontsize=15,fontweight="bold")
        fig.text(.50,-.015,"21 × 21 grid · ideal one-qubit gates · fixed clean calibration · solid: Ising/KT; dashed: PT reference",ha="center",fontsize=9)
        fig.savefig(FIG/f"phase_p{p:.2f}.png",dpi=170,bbox_inches="tight"); plt.close(fig)
    fig,axs=plt.subplots(2,3,figsize=(13,8),layout="constrained")
    for ip,p in enumerate(prep["p"]):
        map_plot(axs[0,ip],prep,labels[ip],f"p = {p:g}")
        im=map_plot(axs[1,ip],prep,np.maximum(f[ip],a[ip]),"Order signal",False)
    fig.colorbar(im,ax=axs[1,:],label="max(F, A)",shrink=.8)
    fig.legend(handles=legend,loc="outside lower center",ncols=3)
    fig.suptitle("A low-noise gate can still erase order in a 254-CNOT circuit",fontsize=17,fontweight="bold")
    fig.savefig(FIG/"preparation_triptych.png",dpi=170,bbox_inches="tight");plt.close(fig)
    # Higher-resolution clean map and a separate fidelity diagnostic.
    fig,axs=plt.subplots(1,3,figsize=(15,4.4),layout="constrained")
    map_plot(axs[0],clean,classify(clean["corr"],t12),"N = 12 · 41 × 41 clean map")
    chi=np.log10(1+clean["chi"])
    im=axs[1].pcolormesh(clean["kappa"],clean["h"]-.025,chi,cmap="viridis",shading="nearest")
    overlay(axs[1]);axs[1].set_title("Independent fidelity response")
    fig.colorbar(im,ax=axs[1],label="log₁₀(1 + χF)",shrink=.8)
    sq=structure_for_plot=order(clean["corr"])[2]
    q=np.argmax(sq,axis=-1)*2/12
    im=axs[2].pcolormesh(clean["kappa"],clean["h"],q,cmap="plasma",vmin=0,vmax=.5,shading="nearest")
    overlay(axs[2]);axs[2].set_title("Dominant correlation momentum")
    fig.colorbar(im,ax=axs[2],label="q*/π; modulation ≠ floating proof",shrink=.8)
    fig.savefig(FIG/"clean_N12_diagnostics.png",dpi=170,bbox_inches="tight");plt.close(fig)
    # Size study: each size uses its own exact Ising calibration, stated openly.
    cuts=json.loads((ROOT/"data/line_cuts.json").read_text())
    rows=[];fig,axs=plt.subplots(1,3,figsize=(13,4),layout="constrained")
    for r in cuts:
        n=r["n"];k=r["kappa"];hh=np.array(r["h"]);cc=np.array(r["corr"])
        ff,aa,_=order(cc);tt=threshold(n);chi=np.array(r["chi"],dtype=float)
        row={"N":n,"kappa":k,"order_boundary":crossing(hh,ff if k<.5 else aa,tt),
             "fidelity_peak":float(hh[np.nanargmax(chi)]-.01)}
        rows.append(row)
        for ax,wanted in zip(axs,[0,.3,.8]):
            if k==wanted:
                ax.plot(hh-.01,chi/n,label=f"N = {n}")
                ax.set(xlabel="Transverse field h",ylabel="χF / N",title=f"κ = {k:g}",xlim=(0,1.4))
                ax.legend(fontsize=9)
    fig.savefig(FIG/"finite_size.png",dpi=170,bbox_inches="tight");plt.close(fig)
    summary["finite_size_boundaries"]=rows
    # ZNE is numerical noise scaling, not a hardware fold claim.
    fig,axs=plt.subplots(1,3,figsize=(13,4.5),layout="constrained")
    map_plot(axs[0],prep,labels[0],"Clean reference")
    map_plot(axs[1],prep,zlabels[0],"Linear ZNE: p=.01 and .02")
    map_plot(axs[2],prep,zlabels[1],"Linear ZNE: p=.05 and .10")
    fig.suptitle("Mitigation recovers part of the signal; high-noise bias remains",fontweight="bold")
    fig.savefig(FIG/"mitigation.png",dpi=170,bbox_inches="tight");plt.close(fig)
    # Boundary shifts: use interpolation only when a downward crossing exists.
    bounds=[]
    for k in [0,.3,.6,.8,1.]:
        j=int(np.argmin(abs(prep["kappa"]-k)))
        for ip,p in enumerate(prep["p"]):
            values=f[ip,:,j] if k<.5 else a[ip,:,j]
            bounds.append({"kappa":k,"p":float(p),"boundary":crossing(prep["h"],values,t8)})
    summary["preparation_boundaries"]=bounds
    # Optional routed readout test. Direct Z measurement uses zero CNOTs.
    cr=clean["corr"];rf,ra,_=order(cr)
    fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout="constrained")
    for ax,k in zip(axs,[.3,.8]):
        j=int(np.argmin(abs(clean["kappa"]-k)))
        for p in [0,.01,.05]:
            ff,aa,_=order(cr[:,j]*attenuation(12,p));y=ff if k<.5 else aa
            ax.plot(clean["h"],y,label=f"Routed parity p={p:g}")
        ax.axhline(t12,color="k",ls="--",lw=.8,label="Fixed threshold")
        ax.set(xlabel="Transverse field h",ylabel="Order signal",title=f"Optional readout test · κ={k:g}",xlim=(0,1.1))
        ax.legend(fontsize=8)
    fig.suptitle("Direct Z readout avoids every routing gate in this test",fontweight="bold")
    fig.savefig(FIG/"readout_control.png",dpi=170,bbox_inches="tight");plt.close(fig)
    controls=json.loads((ROOT/"variational_control/comparison.json").read_text())["rows"]
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout="constrained")
    for ax,name,label in zip(axs,["ferromagnetic","antiphase","paramagnetic"],
                              ["Ferromagnetic point","Antiphase point","Paramagnetic point"]):
        rows=[x for x in controls if x["phase_point"]==name]
        for key,lab,color in [("shallow","28 CNOTs · VQE",COLORS[0]),
                              ("full_amplitude","254 CNOTs · exact loading",COLORS[1])]:
            ax.plot([x["p"] for x in rows],[x[key]["fidelity_to_exact"] for x in rows],
                    "-o",label=lab,color=color,lw=2)
        ax.set(xlabel="Error probability per CNOT p",ylabel="Fidelity to exact ground state",
               ylim=(-.02,1.04),title=label)
        ax.legend(fontsize=8)
    fig.suptitle("Fewer gates preserve the state: same Hamiltonian, same noise",fontweight="bold",fontsize=14)
    fig.savefig(FIG/"circuit_comparison.png",dpi=180,bbox_inches="tight");plt.close(fig)
    with (ROOT/"data/boundaries.csv").open("w",newline="") as fp:
        w=csv.DictWriter(fp,fieldnames=["kappa","p","boundary"]);w.writeheader();w.writerows(bounds)
    (ROOT/"data/summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False))
    print(json.dumps(summary,indent=2))


if __name__=="__main__":run()
