import numpy as np
import math
import sys,os
sys.path.append('../../')
import pickle as pckl
from matplotlib import pyplot as plt
#from tools.mkm_creator import *
#from tools_for_analysis.mkm_creator import *
from tools_for_analysis.FED_tools import plot_FED_with_barrier,read_calculated_data
from tools_for_analysis.intermediates_dict import ads_and_electron

from matplotlib.ticker import FormatStrFormatter


#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


#inputfile = 'C2model_h2o.pkl'
show = True
products=['C2H4','CH3CH2OH','CH4','H2']

potentials=[3.1]
pH = [13]

figsize=(8,8)

alldata=read_calculated_data(None,['100'],1,pklfile='../../parsed_data_for_analysis/parsed_data.pkl',indict=ads_and_electron)

for prod in products:
    ads_and_electron[prod+'_g']['G_vs_pot_100'] = [0,ads_and_electron[prod+'_g']['G']]

if not os.path.exists('FEDs'):
    os.mkdir('FEDs')
if 1:
    plot_FED_with_barrier(ads_and_electron,'100',
                #[['2CO','OCCO','OCCOH','CCO','HCCO','HCCOH','CCH','HCCH','H2CCH','C2H4_g'],
                #['2CO','OCCO','OCCOH','CCO','HCCO','H2CCO','OHCCH2','OCHCH3','H3CCH2O','CH3CH2OH_g'],
                #['2CO','OCCO','OCCOH','CCO','HCCO']],
                [['OCCOH','CCO','HCCO','HCCOH','CCH','HCCH','H2CCH','C2H4_g'],
                ['OCCOH','CCO','HCCO','H2CCO','OHCCH2','OCHCH3','H3CCH2O','CH3CH2OH_g'],
                ['OCCOH','CCO','HCCO']],
                potentials=potentials,pH=[7,13],
                ylim=[-8.9,-0.3],
                view=True,
                annotate_intermediates=False,
                normalize_to_IS=False,
#                check_barriers_with_line=True,
                emphasize_barriers={'pH7':[['HCCO','H2CCO'],['HCCO','HCCOH']],'pH13':[['HCCO','H2CCO'],['HCCOH','CCH']]},
                annotate_reaction_conditions=False,
                linestyles=['-','--'],
                colors=['k','peru','b'],
                figsize=[13,7],
                title='',
                add_potential_response_in_x=True,
                outformat='pdf',
                xticks='nh',
                outdir='FEDs')
if 0:
    # Old style with only one pH included
    plot_FED_with_barrier(ads_and_electron,'100',
                [['CO+CO','OCCO','OCCOH','CCO','HCCO','HCCOH','CCH','HCCH','H2CCH','C2H4_g'],
                ['CO+CO','OCCO','OCCOH','CCO','HCCO','H2CCO','OHCCH2','OCHCH3','H3CCH2O','CH3CH2OH_g'],
                ['CO+CO','OCCO','OCCOH','CCO','HCCO']],
                potentials=potentials,pH=[13],
                ylim=[-4.8,0.8],
                view=True,
                annotate_intermediates=False,
                normalize_to_IS=False,
#                check_barriers_with_line=True,
#                emphasize_barriers={'pH7':[['HCCO','H2CCO'],['HCCO','HCCOH']],'pH13':[['HCCO','H2CCO'],['HCCOH','CCH']]},
                linestyles=['-','--'],
                colors=['k','peru','b'],
                figsize=[12,4],
                title='',
                outformat='pdf',
                outdir='FEDs')
potentials=[3.2]
if 1:
    plt.rcParams["figure.figsize"] = figsize
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    #plt.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plot_FED_with_barrier(ads_and_electron,'100',
                [['OCCOH','HOCCOH'],['OCCOH','CCO'],['OCCOH']],
                potentials=potentials,pH=pH,
                ylim=[-1.8,0.5],
                view=False,
                annotate_intermediates=False,
                normalize_to_IS=True,
#                check_barriers_with_line=True,
                colors=['k','g','r'],
                figsize=figsize,
                title='',
                outformat='pdf',
                outdir='FEDs')

if 1:
    # Single elementary steps
    plot_FED_with_barrier(ads_and_electron,'100',
                [['CCO','HCCO'],['CCO','CCHO'],['CCO','CCOH'],['CCO']],
                potentials=potentials,pH=pH,
                ylim=[-0.6,1.6],
                view=False,
                annotate_intermediates=False,
                normalize_to_IS=True,
#                check_barriers_with_line=True,
                colors=['k','r','tab:orange','g'],
                figsize=figsize,
                outformat='pdf',
                outdir='FEDs',
                title='')

    plot_FED_with_barrier(ads_and_electron,'100',
                [['HCCO','HCCOH'],['HCCO','H2CCO'],['HCCO','HCCHO'],['HCCO']],
                potentials=potentials,pH=pH,
                ylim=[-1.1,0.8],
                view=False,
                annotate_intermediates=False,
                normalize_to_IS=True,
#                check_barriers_with_line=True,
                colors=['k','r','tab:orange','g'],
                figsize=figsize,
                outformat='pdf',
                outdir='FEDs',
                title='')

    plot_FED_with_barrier(ads_and_electron,'100',
                [['H2CCO','OHCCH2'],['H2CCO','H3CCO'],['H2CCO','H2CCOH'],['H2CCO']],
                potentials=potentials,pH=pH,
                ylim=[-1.3,0.6],
                view=False,
                annotate_intermediates=False,
                normalize_to_IS=True,
#                check_barriers_with_line=True,
                colors=['k','r','tab:orange','g'],
                figsize=figsize,
                outformat='pdf',
                outdir='FEDs',
                title='')

