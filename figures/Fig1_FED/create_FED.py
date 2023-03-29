import numpy as np
import math
import sys,os
sys.path.append('../../')
import pickle as pckl
from matplotlib import pyplot
#from tools.mkm_creator import *
#from tools_for_analysis.mkm_creator import *
from tools_for_analysis.FED_tools import plot_diagram,read_calculated_data
from intermediates_dict import ads_and_electron



#inputfile = 'C2model_h2o.pkl'
show = True
products=['C2H4','CH3CH2OH','CH4','H2']


pH = 13
SHE_potential=-0.77
show_reaction_steps=False
RHE_potential = SHE_potential+0.059*pH
c2h4color='b'
ethacolor='peru'
bothc='g'


ylabel='$\Delta$G$^\phi$ / eV'

xlabel='N$_\mathrm{H}$ + $\gamma$'
inputfile=None

figsize=(12,9)

alldata=read_calculated_data(inputfile,'100',1,pklfile='../../parsed_data_for_analysis/parsed_data.pkl',indict=ads_and_electron)
plot_diagram(ads_and_electron,'100',
                show_reaction_steps=show_reaction_steps,
                RHE_potential=RHE_potential,pH=pH,add_field=False,
                ylabel=ylabel,xlabel=xlabel,show=show,yrange=[-1.3,1.5],size=(16,6),
                clean_dict_from_dead_ends=True,add_CO_to_C1s=True,add_title=False,
                name='Figure1.pdf',
                highlights={'HCCHO':bothc, 'OHCCH2':bothc,'CO':bothc,'OCCO':bothc,'OCCOH':bothc,'CCO':bothc,'HCCO':bothc,
                    'OCHCH3':ethacolor,'CH3CH2OH_g':ethacolor,'H3CCH2O':ethacolor,
                    'H2CCH2O':{'self':c2h4color,'O':c2h4color},'O':c2h4color,'OH':c2h4color,'C2H4_g':c2h4color},
                add_potential_response_in_x=True,state_width=0.25)

