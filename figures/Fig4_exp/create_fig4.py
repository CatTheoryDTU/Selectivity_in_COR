#from plot_paper_data import load_pkl_data, _compute_partial_jecsa
import matplotlib.patheffects as path_effects
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from ase import units
import numpy as np
import string,sys
sys.path.append('../../tools_for_analysis/')
from experimental_database_tools import perform_Tafel_analysis,plot_selectivity


plt.rcParams["figure.figsize"] = (8,6.5)
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
markersize=10

Ueqs={'Ethylene':0.08,
      'Ethanol':0.09,
      'Acetate':0.11,
      'n-propanol':0.1,
      'Ethylene+Ethanol+Acetate+n-propanol':0.08}

CRED = '\033[91m'
CEND = '\033[0m'
path_to_db='../../experimental_database/'

finalpaths = [
        "Bertheussen_COR", "Bertheussen_COR-pcCu", "Huang_CO2R",
        "Kanan_COR-ODCu", "Kuhl_CO2", "Wang_COR", "Wang_COR-Cuflower", "Raciti_COR",
        "Jouny_COR", "Luc_COR", "Ma_CO2R",
        "Gregorio_CO2R",
        'Hahn_CO2R','Kanan_COR-ODCuII','Luc_COR_II',
        "Zuettel_CO2R", "Kanan_COR-GDE",
        "Ma_GDE_prelim",
        'Manthiram_CO2R',
        'Varela_CO2R','Yeo_CO2R',#'deLuna_CO2R',
       # "Kanan_CO2R-ODCu", # Taken out because it is quite off the trend of other (too low overpotential)
        ]

folders=[path_to_db+i for i in finalpaths]

emphasized_folders=[]
products = ['Ethylene','Ethanol','Acetate']

general_kw={'emphasized_folders':emphasized_folders,
        'plot_individual_fits':False,
        'mode':None,
        'potscale':'SHE',
        'individual_Tafel_slopes_in_label':False,
        'max_standard_deviation':0.95,
        'fit_CO_and_CO2_separately':False,
        'xlim':[-1.69,-0.95],
#        'minpH':0,
        'color_by_pH':False,
        'show_rsquared': True
        }

fig,ax=plt.subplots(1,1,figsize=(10,7))
studies,folderlabels,all_folders,fits=perform_Tafel_analysis(folders,
         ['+'.join(products)],
        **general_kw,
        verbose=True,
        ax=ax)

plt.legend(bbox_to_anchor=(1.,1),fontsize=8.0,loc='upper left')
plt.close()

bakk_kw=general_kw['mode']
general_kw['mode'] = 'COR'
d1,d2,d3,fitsCO = perform_Tafel_analysis(folders,
        products,#minpH=12,
        **general_kw,
        studies=studies,folderlabels=folderlabels)
general_kw['mode'] = 'CO2R'
d1,d2,d3,fitsCO2 = perform_Tafel_analysis(folders,
        products,
        **general_kw,
        studies=studies,folderlabels=folderlabels)
plt.close()
general_kw['mode']=bakk_kw


fig=plt.figure(figsize=(15,10))

axes=[]
axh=0.2
lowb=0.54
highb=0.76
offset=0.08
axw=0.40

if len(products) == 3:
    axh=0.29
    axes.append(fig.add_axes([offset,offset+2*axh,axw,axh]))
    axes.append(fig.add_axes([offset,offset+axh,axw,axh]))
    axes.append(fig.add_axes([offset,offset,axw,axh]))
    axSCO2= fig.add_axes([2*offset+axw,offset+0.43,axw,0.42])
    axSCO = fig.add_axes([2*offset+axw,offset,axw,0.42])
    axes[0].set_xticks([])
    axes[1].set_xticks([])
    for i in axes:
        i.set_ylabel('log (j$_{\mathrm{ECSA}}$ / mA/cm$^2$)',fontsize=18)

elif len(products) == 2:
    axh=0.29
    axes.append(fig.add_axes([offset,offset+2*axh,axw,axh]))
    axes.append(fig.add_axes([offset,offset+axh,axw,axh]))
    axes.append(fig.add_axes([offset,offset,axw,axh]))
    axSCO2= fig.add_axes([2*offset+axw,offset+0.45,axw,0.42])
    axSCO = fig.add_axes([2*offset+axw,offset,axw,0.42])
    axes[0].set_xticks([])
    axes[1].set_xticks([])
    for i in axes:
        i.set_ylabel('log (j$_{\mathrm{ECSA}}$ / mA/cm$^2$)',fontsize=18)
else:
    axes.append(fig.add_axes([0.22+offset,lowb+axh,axw,axh]))
    axes.append(fig.add_axes([offset,lowb,axw,axh]))
    axes.append(fig.add_axes([axw+offset,lowb+axh,axw,axh]))
    axes.append(fig.add_axes([axw+offset,lowb,axw,axh]))

S = plot_selectivity(fitsCO,colors=['b','peru','r','k'],products=['Ethylene','Ethanol+Acetate'],plot_FE=False,ax=axSCO,mode='COR')
S = plot_selectivity(fitsCO2,colors=['b','peru','r','k'],products=['Ethylene','Ethanol+Acetate'],plot_FE=False,ax=axSCO2,mode='CO2R')


if len(products) == 3:
    axSCO2.set_xlabel('')
    axSCO2.set_xticks([])

else:
    axSCO.set_ylabel('')
    axSCO.set_yticks([])

axSCO2.annotate('(d)',(-1.3,0.97), fontsize=24,color='k',va='top').draggable()
axSCO2.annotate('CO$_2$R',(-1.25,0.97), fontsize=24,fontweight='bold',color='k',va='top').draggable()
axSCO.annotate('(e)',(-1.45,0.975), fontsize=24,color='k',va='top').draggable()
axSCO.annotate('COR',(-1.4,0.97), fontsize=24,fontweight='bold',color='k',va='top').draggable()

for iprod,product in enumerate(products):
      perform_Tafel_analysis(folders,
            product,
            **general_kw,
            ax=axes[iprod],
            studies=studies,folderlabels=folderlabels)
      tag='abcd'[iprod]
      axes[iprod].annotate(f'({tag})',(-1.68,-4.4),va='center',ha='left',fontsize=20)

#Add theoretical selectivities
theosel=np.loadtxt('Theoretical_selectivities.list')

axSCO2.plot(theosel[:,0],theosel[:,1],':b')
axSCO2.plot(theosel[:,0],theosel[:,2],':',color='peru')
axSCO.plot(theosel[:,0],theosel[:,3],':b')
axSCO.plot(theosel[:,0],theosel[:,4],':',color='peru')

legpos=[-1.68,0.55]
legfontsize=20
axSCO2.annotate('Experiment',legpos,arrowprops=dict(arrowstyle="-",relpos=(0,0.6),linestyle='-',linewidth=2),
        xytext=(legpos[0]+0.05,legpos[1]-0.0225),fontsize=legfontsize)
legpos=[legpos[0],legpos[1]-0.1]
axSCO2.annotate('Theory',legpos,arrowprops=dict(arrowstyle="-",relpos=(0,0.6),linestyle=':',linewidth=2),
        xytext=(legpos[0]+0.05,legpos[1]-0.0225),fontsize=legfontsize)

legpos=[-1.18,0.55]
legfontsize=20
axSCO.annotate('Experiment',legpos,arrowprops=dict(arrowstyle="-",relpos=(0,0.6),linestyle='-',linewidth=2),
        xytext=(legpos[0]+0.05,legpos[1]-0.0225),fontsize=legfontsize)
legpos=[legpos[0],legpos[1]-0.1]
axSCO.annotate('Theory',legpos,arrowprops=dict(arrowstyle="-",relpos=(0,0.6),linestyle=':',linewidth=2),
        xytext=(legpos[0]+0.05,legpos[1]-0.0225),fontsize=legfontsize)

axSCO.set_ylim([0.0001,1])
axSCO2.set_ylim([0.0001,1])

plt.savefig('Fig4.pdf')

plt.show()

