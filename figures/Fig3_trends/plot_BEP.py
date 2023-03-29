from ase.io import read,write
from scipy.optimize import curve_fit
import matplotlib
import pickle
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import os,sys
import numpy as np
from general_tools import get_reference_energies,lin_fun,quad_fun
sys.path.append('../../')
from tools_for_analysis.sjm_analyse_tools import plot_BEP_relation,plot_betas,plot_intrinsic_barrier_energies,get_intcolors
from ase import units

#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
markersize=10

home=os.getcwd()

alldata = pickle.load(open('../../parsed_data_for_analysis/parsed_data.pkl','rb'))

fig=plt.figure(figsize=(12,10))

ax_BEP=fig.add_axes([0.08,0.55,0.4,0.42])
ax_beta=fig.add_axes([0.55,0.55,0.4,0.42])
ax_int=fig.add_axes([0.08,0.02,0.4,0.45])

# BEP plot
ax,markers=plot_BEP_relation(alldata,'100',output_potential=3.63,pH=13,
            outdir='BEP/',
            xlabel='$\Delta$E$_{0\mathrm{V}_{\mathrm{RHE}}, \mathrm{pH}=13}$ / eV',
            ylabel='$\Delta$E$^\u2021_{0\mathrm{V}_{\mathrm{RHE}}, \mathrm{pH}=13}$ / eV',
            #ylabel='$\Delta$E$^\u2021$(0V$_{\mathrm{RHE}}$, pH=13)',
            proton_donors=['base'],return_plt=True,figax=(fig,ax_BEP),
            plot_legend=False)

plot_betas(alldata,exclude=['COH-OCCOH','CHO-OCCHO','COCO-OCCO','clean-H','CC-H','H-CCH'],outdir='.',ylabel=r'$\alpha$ / (eV/V)',
        return_plt=True,ax=ax_beta,add_legend=False,markers=markers,xticks=[],plot_meanval_by_type=True)
plot_intrinsic_barrier_energies(alldata,'100',outdir='.',proton_donors=['base'],return_plt=True,figax=(fig,ax_int),ylabel='$\Delta$E$^\u2021_{int}$ / eV',
        markers=markers,xticks=[])

ax_BEP.annotate('(a)',(-1.0,2.2),fontsize=28)
ax_beta.annotate('(b)',(-1.0,0.98),fontsize=28,va='top')

ax_int.annotate('(c)',(-1,1.76),fontsize=28,va='top')


ax_leg=fig.add_axes([0.55,0.02,0.4,0.45],)
fig.patch.set_visible(False)
ax_leg.axis('off')
from matplotlib.pyplot import legend,gca
for ads in markers:
    for toads in markers[ads]:
     if get_intcolors(ads,toads) not in ['k']:
        ax_leg.plot(np.nan,np.nan,markers[ads][toads],color=get_intcolors(ads,toads),label=get_intcolors(ads,toads,return_name=True),markeredgecolor='k',markersize=12)
ax_leg.legend(ncol=2,loc='lower left',fontsize=14,facecolor='none',edgecolor='none',bbox_to_anchor=(0.05,0.0))

ax_leg.set_xlim([0,1])
ax_leg.set_ylim([0,1])
legy=0.880
legx=0.020
legspac=0.5
legyspac=0.07
fontsize=17
ax_leg.plot(legx,legy+legyspac,'s',color='brown',markersize=12)#, markeredgecolor='w')
ax_leg.annotate('C-protonation',(legx+0.05,legy+legyspac),ha='left',va='center',fontsize=fontsize)
ax_leg.plot(legx+legspac,legy+legyspac,'sb',markersize=12)
ax_leg.annotate('O-protonation',(legx+legspac+0.05,legy+legyspac),ha='left',va='center',fontsize=fontsize)
ax_leg.plot(legx,legy,'sy',markersize=12)
ax_leg.annotate('OH-desorption',(legx+0.05,legy),ha='left',va='center',fontsize=fontsize)
ax_leg.plot(legx+legspac,legy,'sg',markersize=12)
ax_leg.annotate('Anionic desorption',(legx+legspac+0.05,legy),ha='left',va='center',fontsize=fontsize)
#plt.tight_layout()
plt.show()
plt.close()
exit()
