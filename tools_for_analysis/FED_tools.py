import numpy as np
import math
import sys,os
import pickle as pckl
from matplotlib import pyplot as plt
from general_tools import get_reference_vibrational_contribution
from ase  import units
#from mkm_creator import *

plt.rcParams["figure.figsize"] = (10,5)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (10,5)
markersize=10


def _clean_dict_from_dead_ends(ads_and_electron,facet):
    allprod=[]

    for ads in reversed(list(ads_and_electron.keys())):
        safe=0
        if 'leads_to' not in ads_and_electron[ads].keys():
            del ads_and_electron[ads]
            continue
        for toads in ads_and_electron[ads]['leads_to']:
            if toads in ads_and_electron.keys():
                if 'E_%s'%facet in ads_and_electron[toads].keys() or toads[-2:] == '_g':
                    safe=1
                    continue
            #if toads in products:
            #        safe =1


        if ads[-2:] == '_g': #in ['H2_g']+products:
            safe=1
        if not safe:
            print(ads+' doesnt seem to lead anywhere')
            del ads_and_electron[ads]
            continue

    for ads in ads_and_electron.keys():
        allprod.extend(ads_and_electron[ads]['leads_to'])
    allprod.append('CO')
    allprod = np.unique(allprod)

    dellist=[]
    for ads in ads_and_electron.keys():
        #if ads[-2:] != '_g' or ads not in allprod:# in allprod and ads not in ['CO_g','H2_g']:
        if ads not in allprod:# in allprod and ads not in ['CO_g','H2_g']:
            print(ads+' doesnt seem to be created')
            dellist.append(ads)
            continue
    for delet in dellist:
        del ads_and_electron[delet]

def add_vibrational_contribution(ads_and_electron,facet,products=[],barrier=False):
    from ase.units import invcm
    from ase.thermochemistry import HarmonicThermo,IdealGasThermo
    from ase.io import read

    for ads in ads_and_electron.keys():
        #If vibrations have been added already continue
        if 'G_vs_pot_%s'%facet in ads_and_electron[ads]: continue
        #If the adsorbate doesn't have an energy on this facet continue
        if 'E_vs_pot_%s'%facet not in ads_and_electron[ads]: continue
        #If no vibrational info is given continue
        if ('vibs_%s'%facet not in ads_and_electron[ads].keys() and
            'free_en_corr_%s'%facet not in ads_and_electron[ads].keys() and
            ads != 'COCO'):
            print('No vibrational information found for %s on %s'%(ads,facet))
            continue

        if ads[-2:] != '_g':
            ads_short_for_vibs=ads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')
            # The vibrations of 2 adjacent *CO is approximated by 2 independent *CO
            # Not ideal but never really used
            if ads == 'COCO':
                #print(ads_and_electron['CO'])
                if 'vibs_%s'%facet in ads_and_electron['CO']:
                  ads_and_electron[ads]['vibs_%s'%facet]=ads_and_electron['CO']['vibs_%s'%facet]+\
                                                       ads_and_electron['CO']['vibs_%s'%facet]
                elif 'free_en_corr_%s'%facet in ads_and_electron['CO']:
                    ads_and_electron[ads]['free_en_corr_%s'%facet] = \
                            ads_and_electron[ads]['free_en_corr_%s'%facet].copy()

            #If the pkl file is used as input
            if 'free_en_corr_%s'%facet in ads_and_electron[ads].keys():
                eng_corr=ads_and_electron[ads]['free_en_corr_%s'%facet]
            else:
                #If catmap input file is used as input
                print(ads)
                df
                gibbs = HarmonicThermo(vib_energies = np.array(ads_and_electron[ads]['vibs_%s'%facet])*invcm)
                eng_corr = gibbs.get_helmholtz_energy(298, verbose=False)
                gas_G_reference=get_reference_vibrational_contribution(ads_short_for_vibs)
                eng_corr-=gas_G_reference
                ads_and_electron[ads]['free_en_corr_%s'%facet] = eng_corr

            #if isinstance(ads_and_electron[ads]['E_%s'%facet],float):
            #    ads_and_electron[ads]['E_%s'%facet] += eng_corr
            #elif isinstance(ads_and_electron[ads]['E_%s'%facet],dict):
            #    for pot in ads_and_electron[ads]['E_%s'%facet].keys():
            #        ads_and_electron[ads]['E_%s'%facet][pot]+=eng_corr

                if ads_and_electron[ads]['E_vs_pot_%s'%facet] is not None:
                    ads_and_electron[ads]['G_vs_pot_%s'%facet] = ads_and_electron[ads]['E_vs_pot_%s'%facet].copy()
                    ads_and_electron[ads]['G_vs_pot_%s'%facet][1] += eng_corr

                else:
                    print(ads+'has no beta!')
        else:
            if 'free_en_corr' in ads_and_electron[ads].keys() and 'G' not in ads_and_electron[ads]:
                ads_and_electron[ads]['G']= ads_and_electron[ads]['E']+ads_and_electron[ads]['free_en_corr']
                continue

            lin = 'nonlinear'
            if ads in ['H2_g','CO_g','CO2_g']:
                lin = 'linear'
                symmetry=1
                pressure=101325.
            elif ads == 'H2O_g':
                pressure=0.035*101325.
                symmetry=2
            else:
                pressure=1
            if ads == 'C2H4_g': symmetry=4
            elif ads == 'CH3CH2OH_g': symmetry=1
            elif ads == 'CH4_g': symmetry=12
            gibbs = IdealGasThermo(vib_energies = np.array(ads_and_electron[ads]['vibs_%s'%facet])*invcm,
                                    geometry=lin,spin=0,symmetrynumber=symmetry,
                                    atoms=read('/Users/geokast/SelectCO2/endstates/gas_geometries/'+ads.rstrip('_g')+'.traj'))

            eng_corr = gibbs.get_gibbs_energy(298.15, pressure=pressure,verbose=False)
            ads_and_electron[ads]['free_en_corr']=eng_corr

            #Change the gas name in order to subtract the reference gas free energies
            ads_for_vibs=ads.replace('H3','HHH').replace('H2','HH').replace('O2','OO').replace('H4','HHHH').rstrip('_g')
            ads_for_vibs=ads_for_vibs.replace('C2','CC')

            #Write Free energy in the dict
            ads_and_electron[ads]['G'] = ads_and_electron[ads]['E'] + eng_corr
            ads_and_electron[ads]['G'] -= get_reference_vibrational_contribution(ads_for_vibs)

        #Add vibrations to transition states
        if ('E_ddag_vs_pot_%s'%facet not in ads_and_electron[ads].keys() or
            ('vibs_ddag_%s'%facet not in ads_and_electron[ads].keys() and
             'free_en_corr_ddag_%s'%facet not in ads_and_electron[ads].keys())):
                continue

        if 'G_ddag_vs_pot_%s'%facet in ads_and_electron[ads].keys():
            #for toads in ads_and_electron[ads]['free_en_corr_ddag_%s'%facet].keys():
            #    #ads_and_electron[ads]['E_ddag_%s'%facet][toads]+=ads_and_electron[ads]['free_en_corr_ddag_%s'%facet][toads]
            #    ads_and_electron[ads]['E_ddag_vs_pot_%s'%facet][toads][1] +=\
            #        ads_and_electron[ads]['free_en_corr_ddag_%s'%facet][toads]
            continue
        for toads in ads_and_electron[ads]['vibs_ddag_%s'%facet].keys():
                gibbs = HarmonicThermo(vib_energies = np.array(ads_and_electron[ads]['vibs_ddag_%s'%facet][toads])*invcm)
                eng_corr = gibbs.get_helmholtz_energy(298, verbose=False)
                if 'G_ddag_vs_pot_%s'%facet not in ads_and_electron[ads].keys():
                    ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet] = {}
                ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet][toads] = ads_and_electron[ads]['E_ddag_vs_pot_%s'%facet][toads].copy()
                ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet][toads][1] += eng_corr
                #ads_and_electron[ads]['E_ddag_%s'%facet][toads] += eng_corr
                #ads_and_electron[ads]['E_ddag_vs_pot_%s'%facet][toads][1] += eng_corr
                toads_short2=toads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')
                ads_short_for_vibs=ads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')

                for letter in ads_short_for_vibs:
                    toads_short2=toads_short2.replace(letter,'',1)
                toads_for_vibs=toads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')

                if toads_short2 == 'H':
                    #TODO: Don't think adding this H2O is correct
                    #ads_short_for_vibs+='HHO'
                    ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet][toads][1]-=0.6325624956374631
                    #ads_and_electron[ads]['E_ddag_%s'%facet][toads]-=0.6325624956374631

                if ads_short_for_vibs == 'CO' and toads == 'OCCO':
                    #For CO to OCCO, CO reference is added twice
                    ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet][toads][1] -=\
                        get_reference_vibrational_contribution(ads_short_for_vibs)

                ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet][toads][1] -=\
                    get_reference_vibrational_contribution(ads_short_for_vibs)
                #ads_and_electron[ads]['E_ddag_%s'%facet][toads] -=\
                #    get_reference_vibrational_contribution(ads_short_for_vibs)


    #TODO: Check if this is consistent and correct! I think this is a relix from before where I did not reference the gas pahsespecies
    if 0:
        for ads in ads_and_electron.keys():
            if 'nHe' in ads_and_electron[ads].keys():
                if 'E_vs_pot_%s'%facet in ads_and_electron[ads].keys():
                    ads_and_electron[ads]['E_vs_pot_%s'%facet][1] -=\
                        ads_and_electron[ads]['nHe']*ads_and_electron['H2_g']['G']
                #if isinstance(ads_and_electron[ads]['E_%s'%facet],float):
                #    ads_and_electron[ads]['E_%s'%facet] -=\
                #        ads_and_electron[ads]['nHe']*ads_and_electron['H2_g']['E']
                #elif isinstance(ads_and_electron[ads]['E_%s'%facet],dict):
                #    for pot in ads_and_electron[ads]['E_%s'%facet].keys():
                #        ads_and_electron[ads]['E_%s'%facet][pot]-=\
                #        ads_and_electron[ads]['nHe']*ads_and_electron['H2_g']['E']




def collect_field_data(alldata,facet,field_data,products=[]):
    thisfacetdata=field_data[facet].copy()
    unused_field_data=thisfacetdata.copy()
    for ads_orig in thisfacetdata.keys():
     ads=ads_orig.lstrip('bdo-').lstrip('md-')
     if ads in alldata.keys() or ads[::-1] in alldata.keys():
        charges=[]
        for charge in thisfacetdata[ads_orig].keys():
            if isinstance(charge,float):
                if 'Erel' in thisfacetdata[ads_orig][charge].keys():
                    charges.append(charge)
                else:
                    print(ads+' %s doesnt work'%charge)

        lowest_charge,highest_charge=min(charges),max(charges)
        #XXX: Should be a better fit
        dedq = (thisfacetdata[ads_orig][lowest_charge]['Erel'] - thisfacetdata[ads_orig][highest_charge]['Erel'])/\
                (thisfacetdata[ads_orig][lowest_charge]['qA'] - thisfacetdata[ads_orig][highest_charge]['qA'])
        if ads not in alldata.keys():
            ads=ads[::-1]
        if 'dedq' not in alldata[ads].keys():
                alldata[ads]['dedq']={}
        #if facet not in alldata[ads]['dedq'].keys():
        #        alldata[ads]['dedq'][facet]=0
        alldata[ads]['dedq'][facet] = dedq
        del unused_field_data[ads_orig]

    else:
        print(ads+' seems to be in the field data but not in catmap input')

    print('Unused field data:',unused_field_data.keys())
    for ads in alldata.keys():
            if 'dedq' not in alldata[ads].keys():
                alldata[ads]['dedq']={}
            if facet not in alldata[ads]['dedq'].keys():
                alldata[ads]['dedq'][facet]=0
                if ads not  in products:
                    print(ads+' doesnt seem to have field response data on facet '+facet)
    #    print(facet, ads, dedq)
    #das
def draw_paths(ads_and_electron,facet,show_reaction_steps=True,RHE_potential=0,pH=0,energy_type='G',
        add_CO_to_C1s=True,add_CO_to_H=True,highlights={},add_potential_response_in_x=False,state_width=0.3):
    enstring='%s_vs_RHE_%s'%(energy_type,facet)

    for ads in ads_and_electron.keys():
        if enstring in ads_and_electron[ads].keys():
         if 'leads_to' in ads_and_electron[ads].keys():
            for toads in ads_and_electron[ads]['leads_to']:
                if toads not in ads_and_electron.keys():
                    #print(toads+' seems not calculated check')
                    continue

                if np.any([i in ads_and_electron[toads].keys()
                          for i in [enstring,energy_type]]):
                    x = ads_and_electron[ads]['nHe']+state_width
                    dx = (1-2*state_width)+(ads_and_electron[toads]['nHe']-ads_and_electron[ads]['nHe']-1)

                    if add_potential_response_in_x:
                        if '%s_vs_pot_%s'%(energy_type,facet) in ads_and_electron[ads]:
                            x+=ads_and_electron[ads]['%s_vs_pot_%s'%(energy_type,facet)][0]
                            dx-=ads_and_electron[ads]['%s_vs_pot_%s'%(energy_type,facet)][0]
                        if '%s_vs_pot_%s'%(energy_type,facet) in ads_and_electron[toads]:
                            dx+=ads_and_electron[toads]['%s_vs_pot_%s'%(energy_type,facet)][0]#-ads_and_electron[ads]['%s_vs_pot_%s'%(energy_type,facet)][0]

                    y = ads_and_electron[ads][enstring][0]*RHE_potential+\
                            ads_and_electron[ads][enstring][1]

                    if add_CO_to_C1s and ads  == 'CO':
                            y *= 2#*ads_and_electron[ads][enstring]

                    if toads[-2:] == '_g':
                        yto=ads_and_electron[toads][energy_type+'_vs_RHE'][0]*RHE_potential+\
                                ads_and_electron[toads][energy_type+'_vs_RHE'][1]
                    else:
                        yto = ads_and_electron[toads][enstring][0]*RHE_potential+\
                            ads_and_electron[toads][enstring][1]
                    dy = yto-y

                    if ads_and_electron[ads]['nHe'] == ads_and_electron[toads]['nHe']:
                        if not add_potential_response_in_x or dx < 0:
                            x = ads_and_electron[ads]['nHe']-state_width+0.05
                            dx=0
#                        elif dx < 0:
#                    if ads_and_electron[toads]['nHe'] - ads_and_electron[ads]['nHe'] > 1:

                    arrowcol='k'
                    arrowwidth=0.25
                    if ads in highlights and toads in highlights:
                        if isinstance(highlights[ads],str):
                            arrowcol=highlights[toads] if isinstance(highlights[toads],str) else highlights[toads]['self']
                            arrowwidth=1.5
                        elif isinstance(highlights[ads],dict):
                            #print(highlights[ads])
                            #dd
                            if toads in highlights[ads]:
                                arrowcol=highlights[ads][toads]
                                arrowwidth=1.5
                            #    print(highlights[ads][toads])
                            #    dd
                        else:
                            dd
                            print(highlights)
                            #dd

                    if toads in ['C2H4_g','CH3CH2OH_g'] and ads_and_electron[toads]['nHe']-ads_and_electron[ads]['nHe'] > 1:
                        # Arrow if more than one PCET is between ads and toads
                        plt.arrow(x,y,dx,dy,head_width=0.05*(1-0.5*RHE_potential),length_includes_head=True,#linestyle='-',
                                linewidth=0.1)
                    else:
                        #Standard arrow
                        print(ads,toads,arrowcol)
                        plt.arrow(x,y,dx,dy,head_width=0.05*(1-0.5*RHE_potential),length_includes_head=True,color=arrowcol,linewidth=arrowwidth)
                        #plt.arrow(x,y,dx,dy,head_width=0.05*(1-0.5*RHE_potential),length_includes_head=True,linewidth=arrowwidth)


def  add_CHE_and_energy_vs_RHE(ads_and_electron,facet,show_reaction_steps=True,
        add_field=False,PZC=None,Capacitance=None,pH=None,T=298.,V0SHE=4.4,add_CO_to_C1s=True,add_CO_to_H=True):

    #if 'G_vs_pot_%s'%facet not in ads_and_electron['CO'].keys():
    add_vibrational_contribution(ads_and_electron,facet,barrier=True)

    for entype in ['E','G']:
        enstring=entype+'_vs_pot_%s'%facet
        for ads in ads_and_electron.keys():
            if (enstring not in ads_and_electron[ads].keys() and
                ads[-2:] != '_g'): continue
            if 'nHe' not in ads_and_electron[ads].keys(): continue
            nhe = ads_and_electron[ads]['nHe']

            #Gas  phase species:
            if ads[-2:] == '_g':
                print(ads, ads_and_electron[ads].keys())
                if entype in ads_and_electron[ads]:
                    ads_and_electron[ads][entype+'_vs_RHE']=np.array([nhe,ads_and_electron[ads][entype]])
                continue

#            print(ads_and_electron[ads][enstring])
            e_0V = np.array(ads_and_electron[ads][enstring])@[-np.log(10)*units.kB*T*pH+V0SHE,1]
            ads_and_electron[ads][entype+'_vs_RHE_%s'%facet]=\
                np.array([ads_and_electron[ads][enstring][0]+nhe,e_0V])

            #For C1 path an extra CO is added in the energy for consistency with C2 pathway
            #TODO: Should maybe be somewhere else
            if add_CO_to_C1s:
             if  ads in ['CHO','COH','CHOH','CH','CH2','CH3','C','CHOH','CH2OH']:
                #TODO: DO not remember what the following two lines did
                #if show_reaction_steps:
                #    e+=RHE_potential+ads_and_electron['CO'][energy_type+'_vs_RHE']
                #else:
                ads_and_electron[ads][entype+'_vs_RHE_%s'%facet][1]+=ads_and_electron['CO'][entype+'_vs_RHE_%s'%facet][1]
            if add_CO_to_H and ads == 'H':
                ads_and_electron[ads][entype+'_vs_RHE_%s'%facet][1]+=ads_and_electron['CO'][entype+'_vs_RHE_%s'%facet][1]
            #For O and OH ethylene is added for consistency with C2 pathway
            elif ads  in ['O','OH']:
                ads_and_electron[ads][entype+'_vs_RHE_%s'%facet][1]+=ads_and_electron['C2H4_g'][entype]

def plot_diagram(ads_and_electron,facet,show_reaction_steps=True,
        show_boltzmann_weights=False,RHE_potential=0,pH=0,add_field=False,
        ylabel=None,xlabel='H$^+$ and e$^-$ transferred',
        show=False,name='FEDs/Free_energies_%s_%1.2fV.pdf',
        yrange=None,size=(10,5),energy_type='G',clean_dict_from_dead_ends=False,
        add_CO_to_C1s=True,add_CO_to_H=True,add_title=True, add_potential_response_in_x=False,
        highlights={},state_width=0.3,alpha_not_highlighted=0.7):

    if ylabel is None:
        ylabel=r'$\Delta$ %s$^\phi$ / eV'%energy_type

    plt.rcParams["figure.figsize"] = size
    plt.rcParams["xtick.labelsize"] = 25
    plt.rcParams["ytick.labelsize"] = 25
    print('-'*13+'\nPlotting diagram for %s'%facet)
    enstring='%s_vs_RHE_%s'%(energy_type,facet)

    if clean_dict_from_dead_ends:
        _clean_dict_from_dead_ends(ads_and_electron,facet)
    if enstring not in ads_and_electron['CO'].keys():
        add_CHE_and_energy_vs_RHE(ads_and_electron,facet,show_reaction_steps,
                add_field,pH=pH,add_CO_to_C1s=add_CO_to_C1s,add_CO_to_H=add_CO_to_H)

    draw_paths(ads_and_electron,facet,show_reaction_steps,RHE_potential,pH,energy_type,add_CO_to_C1s=add_CO_to_C1s,
            add_CO_to_H=add_CO_to_H,highlights=highlights,add_potential_response_in_x=add_potential_response_in_x,state_width=state_width)

    maxNHE=0
    for ads  in ads_and_electron.keys():
        if (enstring not in ads_and_electron[ads].keys() and
            ads[-2:] != '_g'):
            print('Energy of %s not found'%ads)
        elif 'nHe' in ads_and_electron[ads].keys():
            if ads[-2:] == '_g':
                e=ads_and_electron[ads][energy_type+'_vs_RHE']@[RHE_potential,1]
            else:
                e=ads_and_electron[ads][enstring][0]*RHE_potential+ads_and_electron[ads][enstring][1]

            nhe=ads_and_electron[ads]['nHe']

            if add_potential_response_in_x:
                if energy_type+f'_vs_pot_{facet}' in ads_and_electron[ads]:
                    nhe+=ads_and_electron[ads][energy_type+f'_vs_pot_{facet}'][0]

            if show_boltzmann_weights:
                if nhe not in boltz_sorted_by_nhe.keys():
                    boltz_sorted_by_nhe[nhe]=[]
                boltz_sorted_by_nhe[nhe].append(np.exp(-e/kT))

            color='k'
            alpha=1
            if ads in highlights:
                color=highlights[ads] if isinstance(highlights[ads],str) else highlights[ads]['self']
            elif len(highlights.keys()):
                alpha=alpha_not_highlighted


            if ads not in  ['CO','H2_g','H2O_g','CO_g']:#,'CHO','COH']:
                    output = np.array([[nhe-state_width,e],[nhe+state_width,e]])

                    if 'free_en_corr' in ads_and_electron[ads].keys():
                        if not ads_and_electron[ads]['free_en_corr']:
                            adsstring = ads.rstrip('_g')+'*'
                        else:
                            adsstring = ads.rstrip('_g')
                    elif 'free_en_corr_%s'%facet in ads_and_electron[ads].keys():
                        if not ads_and_electron[ads]['free_en_corr_%s'%facet]:
                            adsstring = ads+'*'
                        else:
                            adsstring = ads
                    elif 'vibs_%s'%facet in ads_and_electron[ads].keys():
                        if len(ads_and_electron[ads]['vibs_%s'%facet]) == 0:
                            adsstring = ads+'*'
                        else:
                            adsstring = ads

                    else:
                            adsstring = ads+'*'

                    #Adding CO to C1 products for consistency
                    if ads in ['CHO','COH','CHOH','CH','CH2','CH3','C','CHOH','CH2OH','CH4_g']:
                        line='--'
                        if add_CO_to_C1s: outstring=adsstring+'+CO'
                        else: outstring=adsstring
                    elif ads == 'H' and add_CO_to_H: outstring,line=adsstring+'+CO','--'
                    #Adding Ethylene to O and OH
                    elif ads in ['O','OH']:outstring,line=adsstring+'+C$_2$H$_{4(\mathrm{g})}$','-'

#                        outstring,line=adsstring,'-'
                    elif any([i.isnumeric() for i in ads[1:]]):
                        outstring=ads[0]
                        for letter in ads[1:]:
                            if letter.isnumeric():
                                outstring+=f'$_{letter}$'
                            else:
                                outstring+=letter
                    else: outstring,line=adsstring,'-'
                    if outstring == 'CH$_3$CH$_2$OH_g':
                        outstring = 'C$_2$H$_5$OH_l'
                    outstring=outstring.replace('_g','$_{(\mathrm{g})}$')
                    outstring=outstring.replace('_l','$_{(\mathrm{l})}$')



                    plt.plot(output[:,0],output[:,1],line,linewidth=3,color=color,alpha=alpha)
                    plt.annotate(outstring,xy=(nhe,e),ha='center',va='bottom',color=color,bbox={'facecolor':'w','edgecolor':'none','pad':-1.5},alpha=alpha).draggable()

            elif ads in ['CO']:
                if add_CO_to_C1s:
                    output = np.array([[nhe-state_width,e*2],[nhe+state_width,e*2]])
                    plt.annotate('2'+ads, xy=(0,e*2+0.02),ha='center',va='bottom',bbox={'facecolor':'w','edgecolor':'none','pad':-1.5},color=color).draggable()
                else:
                    output = np.array([[nhe-state_width,e],[nhe+state_width,e]])
                    plt.annotate(ads, xy=(nhe,e),ha='center',va='bottom',color=color).draggable()
                plt.plot(output[:,0],output[:,1],'-k',linewidth=3,color=color)
            if nhe > maxNHE: maxNHE=nhe

        else:
            print('Seems like adsorbate %s does not have the "nHe" tag'%ads)

    if show_boltzmann_weights:
        boltz_sorted_by_nhe={}
        for nhe in boltz_sorted_by_nhe.keys():
            part_fct=np.sum(boltz_sorted_by_nhe[nhe])
            for ist,state in enumerate(boltz_sorted_by_nhe[nhe]):
                exist_prob=state/part_fct
                try:
                    order_of_mgn=int(math.floor(math.log10(abs(exist_prob))))
                except ValueError:
                    order_of_mgn=13
                leading_value = exist_prob/(10**order_of_mgn)
                en=-np.log(boltz_sorted_by_nhe[nhe][ist])*kT
                if nhe > -0.1:#not in ['CO']:
                    plt.text(nhe-0.5,en,'%ie%i'%(
                        leading_value,order_of_mgn))
                else:
                    plt.text(nhe-0.5,en*2,'%ie%i'%(
                        leading_value,order_of_mgn))


    if not yrange:
        yrange=[-2.5+8*RHE_potential,1.5]

    plt.ylim(yrange)
    plt.xlim([-1.5*state_width,maxNHE+1.7*state_width])

    SHE_potential=RHE_potential-0.059*pH
    if add_title:
        plt.title('Facet: %s, pH=%i, V$_{SHE}$=%1.2f  V$_{RHE}$=%1.2f'%(facet,pH,SHE_potential,RHE_potential))
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    #plt.annotate('Cu({facet})\n{SHE_potential:1.2f}{shepot}\npH {pH:%i}',(0,-0.6),fontsize=22).draggable()
    plt.annotate('Cu(%s)\n%1.2fV$_{\mathrm{SHE}}$\npH%i'%(facet,SHE_potential,pH),
            (-state_width,yrange[0]),va='bottom',ha='left',fontsize=22).draggable()
#    plt.annotate('%1.2fV$_{SHE}$'%(RHE_potential-0.059*pH),(0,-1.0),fontsize=22).draggable()
#    plt.annotate('pH %i'%pH,(0,-1.4),fontsize=22).draggable()
    plt.tight_layout()
    if show:
        plt.show()
    else:
        if add_field:
            plt.savefig(name%(facet,RHE_potential))
        else:
            plt.savefig(name%(facet,RHE_potential))
    plt.close()

def read_calculated_data(inputfile,facets=None,start_from_pkl=False,
                         indict=None,substrates=['Cu'],pklfile='results/parsed_data.pckl',
                         add_field=False,field_data=None,V0SHE=4.4,PZC=None,Capacitance=None):

    if not indict:
            sys.path.append('/Users/geokast/SelectCO2/endstates')
            from intermediates_dict import ads_and_electron

    else:
        ads_and_electron=indict
        #ads_and_electron={}

    if start_from_pkl:
        print('Reading data from %s'%pklfile)
        alldata=pckl.load(open(pklfile,'rb'))
        #if not indict:
        for ads in alldata.keys():
                ads_short=ads.lstrip('md-').lstrip('bdo-')
                if ads_short in ads_and_electron.keys():
                 ads_and_electron[ads_short].update(alldata[ads])
        return

    print('Reading data from %s'%inputfile)

    inlines = open(inputfile,'r').readlines()[1:]
    if isinstance(facets,str):
        facets=[facets]

    #Set up dictionary
    for facet in facets:
      for line in inlines:
         ads = line.split()[2]
         if line.split()[0] in substrates+['None']:
            if line.split()[1] == 'gas':
                ads+='_g'

            if ads in ['CO2_g','CH3COOH_g','CH2CO_g','OH_g']: continue
            if not indict and ads not in ads_and_electron.keys():
                    ads_and_electron[ads] = {}

      if add_field:
        collect_field_data(ads_and_electron,facet,field_data)

    for facet in facets:
      for line in inlines:
         ads = line.split()[2]
         if line.split()[0] not in substrates+['None']: continue
         if line.split()[1] not in [facet,'gas']: continue

         if line.split()[1] == 'gas':
                ads+='_g'

                if ads in ads_and_electron.keys():
                #if ads in ['CO2_g','CH3COOH_g','CH2CO_g','OH_g']: continue
                    ads_and_electron[ads]['E'] = float(line.split()[3])

         #elif line.split()[1] in [facet,'gas']:
         elif ads in ads_and_electron.keys():
             ads_and_electron[ads]['E_%s'%facet] = float(line.split()[3])
             #TODO: Here add field potential dependence!
             if add_field:
                 if not PZC or not Capacitance:
                     print('PZC and Capacitance have to be given to the read function!')
                     ads_and_electron[ads]['E_vs_pot_%s'%facet] = np.array([0,ads_and_electron[ads]['E_%s'%facet]])
                 else:
                     dedphi = ads_and_electron[ads]['dedq'][facet]*Capacitance
                     if PZC < 3: PZC=PZC+V0SHE
                     offset = ads_and_electron[ads]['E_%s'%facet]-PZC*dedphi
                     #print(ads,dedphi)
                     ads_and_electron[ads]['E_vs_pot_%s'%facet] = np.array([dedphi,offset])

             else:
                 ads_and_electron[ads]['E_vs_pot_%s'%facet] = read_beta_from_catmap_input(line)

         if line.split()[1] in [facet,'gas']:
             if ads in ads_and_electron.keys():
                 freq_inline=[None,None]
                 for isplit,splitline in  enumerate(line.split()):
                     if splitline[0] == '[':
                         freq_inline[0]=isplit
                     elif splitline[-1] == ']':
                         freq_inline[1]=isplit+1
                         break

                 if None not in freq_inline:
                     frequencies = [float(vib.replace(',','').replace('[','').replace(']',''))
                             for vib in line.split()[freq_inline[0]:freq_inline[1]]]

                 else:
                     print('No frequencies given for '+ads)
                     frequencies=[]

                 ads_and_electron[ads]['vibs_%s'%facet] = frequencies

             elif '-' in ads:
                 barads=ads.split('-')[0]
                 if barads not in ads_and_electron.keys() and barads != 'OC': continue
                 if ads.split('-')[1] == 'ele':
                     if barads == 'COCO':
                         bartoads='OCCO'

                     elif barads.replace('O','',1).replace('H','',1) in ads_and_electron.keys():
                         bartoads=barads.replace('O','',1).replace('H','',1)
                     elif barads[::-1].replace('O','',1).replace('H','',1) in ads_and_electron.keys():
                         bartoads=barads[::-1].replace('O','',1).replace('H','',1)
                         bartoads=bartoads[::-1]
                     else:
                         print('Could not determine the product of %s'%ads)
                         bartoads=ads

                 elif ads.split('-')[1] == 'H2O':
                     if 'H'+barads in ads_and_electron.keys():
                         bartoads='H'+barads
                     elif barads == 'OC':
                         barads,bartoads='CO','CHO'
                     elif barads == 'CO':
                         barads,bartoads='CO','COH'
                     elif barads+'H' in ads_and_electron.keys():
                         bartoads=barads+'H'
                     else:
                         print('Could not determine the product of %s'%ads,barads)
                         bartoads=ads

                 else:
                     print(ads,' special case check')
                     bartoads=barads[::-1]+barads


                 if 'E_ddag_%s'%facet not in ads_and_electron[barads].keys():
                     ads_and_electron[barads]['E_ddag_%s'%facet]={}
                     ads_and_electron[barads]['E_ddag_vs_pot_%s'%facet]={}
                     ads_and_electron[barads]['vibs_ddag_%s'%facet]={}
                 ads_and_electron[barads]['E_ddag_%s'%facet][bartoads] = float(line.split()[3])
                 ads_and_electron[barads]['E_ddag_vs_pot_%s'%facet][bartoads] =                                read_beta_from_catmap_input(line)

                 freq_inline=[None,None]
                 #da
                 for isplit,splitline in  enumerate(line.split()):
                     if splitline[0] == '[':
                         freq_inline[0]=isplit
                     elif splitline[-1] == ']':
                         freq_inline[1]=isplit+1
                         break

                 if None not in freq_inline:
                     frequencies = [float(vib.replace(',','').replace('[','').replace(']',''))
                             for vib in line.split()[freq_inline[0]:freq_inline[1]]]
                 #ads_and_electron[barads]['vibs_ddag_%s'%facet][bartoads] = read_vibrational_frequencies(ads,   line,vibfile,facet)

                 ads_and_electron[barads]['vibs_ddag_%s'%facet][bartoads] = frequencies
             else:
                 print(ads+' should be added to the dict')
    return ads_and_electron



def read_data_from_catmap_pckl(inputfile):
    all_data=pckl.load(open(inputfile,'rb'))
    #for facet in facets:
    for ads in all_data['species_definitions'].keys():
            if ads.split('_')[-1] == facet:
                    if ads.split('_')[0] in ads_and_electron.keys():
                        ads_and_electron[ads.split('_')[0]]['E'] = all_data['species_definitions'][ads]['formation_energy'][0]

    return ads_and_electron

def facet_difference(ads_and_electron,facets,products=[],ref_facet='100',sort_by_energy=False,specific_adsorbates=[],absolute_energies=False,equalize_adsorbate=None,
        figsize=(6,4)):
#    ads_and_electron=read_data_from_catmap_input(inputfile,facets=facets)
#    print(ads_and_electron)
#    for facet in facets:
#                add_vibrational_contribution(ads_and_electron,facet,products)

    xlabel_list={}
    colors={'100':'b','111':'k','211':'g','110':'r'}
    dE={}
    for facet in facets:
        if facet not in dE.keys():
            dE[facet]=[]
            xlabel_list[facet]=[]
        #if facet == ref_facet:
        #    continue
        for  nHE in range(8):
            for iads,ads in enumerate(ads_and_electron.keys()):
                if (ads in products or
                    'nHe' not in ads_and_electron[ads].keys()):
                    continue

                if ads_and_electron[ads]['nHe'] != nHE:
                    continue

                if not absolute_energies:
                    if 'E_%s'%ref_facet not in ads_and_electron[ads].keys():
                        print(ads+' is not in the reference facet!')
                        continue

                    ref_E=ads_and_electron[ads]['E_%s'%ref_facet]

                if 'E_%s'%facet not in ads_and_electron[ads].keys():
                    print(ads+' is not present on %s!'%facet)
                    continue

                E_rel=ads_and_electron[ads]['E_%s'%facet]
                if not absolute_energies:
                    E_rel-=ref_E

                #if abs(E_rel) < 1e-5:
                #    continue
                dE[facet].append(E_rel)

                if ads in ['O','OH']:
                    xlabel_list[facet].append(r'C$_2$H$_4$ + %s'%ads)
                else:
                    xlabel_list[facet].append(ads)
            dE[facet].append(np.nan)
            xlabel_list[facet].append('')

    fig_all=plt.figure()
    fig_all.dpi = 300
    ax_all=fig_all.add_subplot(1,1,1)
    for ifac,facet in enumerate(facets):
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        if facet == ref_facet:
            continue
        dE[facet]=np.array(dE[facet])

        if sort_by_energy:
            sorti=np.argsort(dE[facet])
            dE[facet]=dE[facet][sorti]
            xlabel_list2=np.array(xlabel_list[facet])[sorti]
        else:
            xlabel_list2=np.array(xlabel_list[facet])
        ax.plot(dE[facet],'o')

        ax.axhline(y=0,linestyle='--',color='k')

        for i,E in enumerate(dE[facet]):
            if  np.isnan(E):
                ax.axvline(x=i,linestyle='--',color='k',linewidth=0.5)

        ax.set_xticks(np.arange(len(xlabel_list2)))
        ax.set_xticklabels(xlabel_list2,rotation='vertical')

        ax.set_title("%s relative to %s"%(facet,ref_facet))

        if absolute_energies:     ax.set_ylabel(r'$\Delta$ G [eV]')
        else:     ax.set_ylabel(r'$\Delta\Delta$ G [eV]')

        fig.tight_layout()
        fig.savefig('results/Binding_energies_%s_rel_%s.pdf'%(facet,ref_facet))

        ax_all.plot(range(len(dE[facet])),dE[facet],'o',color=colors[facet])
        #plt.close()


    ax_all.axhline(y=0,linestyle='--',color='k')
    ax_all.set_xticks(np.arange(len(xlabel_list2)))
    ax_all.set_xticklabels(xlabel_list2,rotation='vertical')
    ax_all.set_title("All facets relative to %s"%(ref_facet))
    for i,E in enumerate(dE[facets[0]]):
        if  np.isnan(E):
                ax_all.axvline(x=i,linestyle='--',color='k',linewidth=0.5)

    if absolute_energies:     ax_all.set_ylabel(r'$\Delta$ G [eV]')
    else:     ax_all.set_ylabel(r'$\Delta\Delta$ G [eV]')

    fig_all.tight_layout()
    fig_all.savefig('results/Binding_energies_all_rel_%s.pdf'%(ref_facet))
    fig_all.savefig('results/Binding_energies_all_rel_%s.png'%(ref_facet),transparent=True)

    for ifac,facet in enumerate(facets):
        if facet == ref_facet:
            continue

    if specific_adsorbates:
        plt.rcParams["figure.figsize"] = figsize
        for ifac,facet in enumerate(facets):
        #    if facet == ref_facet:
        #        continue
            dE[facet]=np.array(dE[facet])
            counter=0
            for iads,ads in enumerate(specific_adsorbates):
                for i,label in enumerate(xlabel_list[facet]):
                    if label == ads:
                        if equalize_adsorbate:
                            dE[facet][i]-=ads_and_electron[equalize_adsorbate]['E_%s'%facet]
                        plt.plot([iads-0.25,iads+0.25],[dE[facet][i],dE[facet][i]],'-',color=colors[facet])

                        if ads == equalize_adsorbate:
                            plt.plot([iads-0.25,iads+0.25],[dE[facet][i],dE[facet][i]],'-',color='k')
            plt.plot(np.nan,np.nan,'-',color=colors[facet],label=facet)

        specific_name='-'.join(specific_adsorbates)
        for iads,adso in enumerate(specific_adsorbates):
            if adso in ['O','OH']:
                specific_adsorbates[iads]=r'C$_2$H$_4$ + %s'%adso

        plt.xticks(np.arange(3),specific_adsorbates,rotation='vertical',fontsize=14)
        if absolute_energies:     plt.ylabel(r'$\Delta$ G [eV]')
        else:     plt.ylabel(r'$\Delta\Delta$ G [eV]')
        plt.legend()
        plt.tight_layout()
    #    plt.show()
        plt.savefig('results/Binding_of_%s.pdf'%(specific_name))
        plt.close()
        plt.rcParams["figure.figsize"] = (10,5)

def plot_FED_with_barrier(alldata,facets,included_steps,potentials=[],pH=[13],
        ylim=[-1.2,0.6],view=False, annotate_intermediates=True,
        energy_type='G',proton_donor='base',V0_SHE=4.4,normalize_to_IS=False,
        figsize=(15,5),title=None,annotate_reaction_conditions=True,
        check_barriers_with_line=False,outdir='results',outformat='png',
        labelsize=None,colors=['k','b','r','g','y'],emphasize_barriers=[],
        linestyles=['-','--','-.',':'],return_plt=False,xticks=None,
        add_potential_response_in_x=False,xlabel='N$_\mathrm{H}$ + $\gamma$',
        include_kinetics=True,
        plotter=None):

    if plotter is not None:
        plt=plotter
    else:
        from matplotlib import pyplot as plt
    plt.rcParams["figure.figsize"] = figsize


    if not potentials:
        print('No potentials where given for FED with barriers')
        return

    if isinstance(included_steps[0],str):
        included_steps=[included_steps]

    if isinstance(facets,str):
        facets=[facets]

    if len(facets) > 1:
        # If more than one facet is included the facets have different colors
        # and the mechanism different linestyles
         colors=list(reversed(colors[:len(facets)]))
         linestyles=list(reversed(linestyles[:len(included_steps)]))
         print('Several facets have been given, they will have varying '
               'colors and the mechanism will change linestyle')

    else:
        # If only one facet is included the mechanism have different colors
        # and the potentials different linestyles
         colors=list(reversed(colors[:len(included_steps)]))
         #linestyles=list(linestyles[0])*len(included_steps)
         linestyles=list(linestyles)*len(included_steps)
         if len(list(linestyles)) >= len(potentials):
             linestyles=list(linestyles)*len(included_steps)
         print('Only one facet has been given the mechanisms will have varying '
               'colors and the potentials will change linestyle')

    _xticks=[]
    all_ens=[]
    lowest_xval,longest_mech=100,0
    for ifac,facet in enumerate(facets):
     #print(len(facets),ifac,colors)
     enstring='%s_vs_pot_%s'%(energy_type,facet)
     barenstring='%s_ddag_vs_pot_%s'%(energy_type,facet)
     for imech,single_mech_steps_in in enumerate(included_steps):
        if len(single_mech_steps_in) > longest_mech: longest_mech=len(single_mech_steps_in)
        #If facets are compared color by facet otherwise by mechanism
        if len(facets) > 1:
            color = colors[ifac]
    #        print(color)
        elif len(included_steps) > 1:
            color = colors[imech]
        else:
            color=colors[0]

        single_mech_steps=single_mech_steps_in
        added_ads=[None]*len(single_mech_steps_in)
        #If C2  is plotted in the same plot as C1 intermediates the name is given
        # with a plus and will be split here
        if np.any([['+' in i for i in single_mech_steps_in]]):
            single_mech_steps = [i.split('+')[0] for i in single_mech_steps_in]
            for istep,step in enumerate(single_mech_steps_in):
                if step.split('+')[-1] != single_mech_steps_in[istep]:
                    added_ads[istep] = step.split('+')[-1]
                else:
                    added_ads.append(None)

        #Recognize if adsorbate is given with a 2 infront e.g. 2CO
        if np.any([[i[0] == '2' for i in single_mech_steps_in if len(i)]]):
            single_mech_steps = [i.lstrip('2') for i in single_mech_steps]
            for istep,step in enumerate(single_mech_steps_in):
                if not len(step): continue
                if step[0] == '2':#single_mech_steps[istep]:
                    added_ads[istep] = step.lstrip('2')#split('+')
                else:
                    added_ads.append(None)

        if enstring not in alldata[single_mech_steps[0]].keys():
            add_CHE_and_energy_vs_RHE(alldata,facet,pH=pH)
            #add_vibrational_contribution(alldata,facet,barrier=[])


        for iph,ph in enumerate(pH):
         if isinstance(proton_donor,dict):
            pdonor=proton_donor[ph]
         else:
             pdonor=proton_donor


         for ipot,potential in enumerate(potentials):
            rhe_pot=potential-(V0_SHE-0.059*ph)
            IS_normalization=0
            for istep,current_intermediate in enumerate(single_mech_steps):
#                current_intermediate=single_mech_steps[istep]
                if istep < len(single_mech_steps)-1:
                    next_intermediate,inext_step=single_mech_steps[istep+1],istep+1
                    if next_intermediate == '':
                        next_intermediate=single_mech_steps[istep+2]
                        inext_step=istep+2
                if current_intermediate not in alldata.keys():
                    print('Could not find the intermediate ',current_intermediate)
                    continue
                if enstring not in alldata[current_intermediate]:
                    print('Intermdiate %s does not seem to have an energy with name %s'%(current_intermediate,enstring))
                    continue

#                print(alldata[current_intermediate].keys(),current_intermediate)
                _xticks.append(int(alldata[current_intermediate]['nHe']))

                En_IS=np.poly1d(alldata[current_intermediate][enstring])(potential)
                #Add CHE to endstate
                En_IS+=alldata[single_mech_steps[istep]]['nHe']*rhe_pot

                if added_ads[istep] is not None:
                    En_IS+=np.poly1d(alldata[added_ads[istep]][enstring])(potential)+\
                            alldata[added_ads[istep]]['nHe']*rhe_pot

                if normalize_to_IS and istep==0:
                    IS_normalization=-En_IS

                if current_intermediate == 'CO' and 'OCCO' in single_mech_steps:
                    if barenstring not in alldata['CO']:
                        alldata['CO'][barenstring] = {}

                    if barenstring in alldata['CO']:
                        if 'OCCO' not in alldata['CO'][barenstring].keys():
                            alldata['CO'][barenstring]['OCCO']= alldata['COCO'][barenstring]['OCCO']



                if annotate_intermediates:
                    stepout=current_intermediate
                    if added_ads[istep] is not None:
                        stepout+='+'+added_ads[istep]
                        if current_intermediate == added_ads[istep]:
                            stepout = '2'+current_intermediate

                    if labelsize is None:
                        plt.annotate(stepout,xy=(inext_step,ylim[1]-0.1*(imech+1)),color=color).draggable()
                    else:
                        plt.annotate(stepout,xy=(inext_step,ylim[1]-0.1*(imech+1)),color=color,fontsize=labelsize).draggable()

                #if istep == len(single_mech_steps)-1:
                #    plt.plot([istep+0.75,istep+1.25],[En_IS+IS_normalization,En_IS+IS_normalization],'-'+colors[imech])
                #    continue
                #print(istep,current_intermediate)
                if istep < len(single_mech_steps)-1:
                    if next_intermediate == '':
                        print('Skipping the step %i'%inext_step)
                    elif next_intermediate not in alldata.keys():
                        print('Could not find the intermediate for FS:',next_intermediate)
                        continue
                    else:
                        if enstring not in alldata[next_intermediate]:
                            print('Intermdiate %s_%s does not seem to have an energy'%(next_intermediate,facet))
                            continue

                    print(current_intermediate,next_intermediate)
                    En_FS=np.poly1d(alldata[next_intermediate][enstring])(potential)+\
                            alldata[next_intermediate]['nHe']*rhe_pot

                    if added_ads[inext_step] is not None:
                        En_FS+=np.poly1d(alldata[added_ads[inext_step]][enstring])(potential)+\
                                        alldata[added_ads[inext_step]]['nHe']*rhe_pot


                    #If no  barrier at all has  been calculated from the current intermediate
                    #Draw a straight line to the next intermediate
                    Eddag=None
                    if barenstring not in alldata[current_intermediate].keys() or not include_kinetics:
                        #XXX: IS the arrow below ever needed?
#                        plt.arrow(istep+1.25,En_IS,0.5,En_FS-En_IS,head_width=0.0,length_includes_head=True, linewidth=0.1,color='r')
                        pass
                    #For multiple H adsorptions
                    elif next_intermediate in ['2H'] and current_intermediate in ['H']:
                        Eddag=np.poly1d(alldata['clean'][barenstring][next_intermediate]['base'])(potential)+alldata['H'][enstring]

                    #For a chemical step
                    elif next_intermediate in ['H','OCCO','CO2']:
                        #Eddag=np.poly1d(alldata['clean'][barenstring]['H']['base'])(potential)
                        Eddag=np.poly1d(alldata[current_intermediate][barenstring][next_intermediate]['base'])(potential)
                    #If the barrier between IS and FS has been calculated
                    elif next_intermediate in alldata[current_intermediate][barenstring]:
                     toads=next_intermediate
                     if pdonor in alldata[current_intermediate][barenstring][toads]:
                        Eddag=np.poly1d(alldata[current_intermediate][barenstring][toads][pdonor])(potential)


                    #print(facet,Eddag,single_mech_steps[istep])

                     #Add CHE to barriers (alkaline has CHE like IS, acid has CHE like FS)
                    if Eddag is not None:
                        Eddag+=alldata[current_intermediate]['nHe']*rhe_pot
                        if added_ads[istep] is not None and added_ads[inext_step] is not None:
                                Eddag+=np.poly1d(alldata[added_ads[istep]][enstring])(potential)


                        if pdonor == 'acid' and next_intermediate not in ['OCCO','CO2']:
                                Eddag+=rhe_pot

                        #Connect the states by lines (TODO: Maybe polynoms?)

                #plot everything
                #Plot IS
                if len(facets) > 1:
                    linestyle=linestyles[imech]
                elif len(potentials) > 1:
                    linestyle=linestyles[ipot]
                elif len(pH):
                    linestyle=linestyles[iph]

                current_xval,next_xval=istep,inext_step
                if add_potential_response_in_x:
                    current_xval=alldata[current_intermediate][enstring][0]+alldata[current_intermediate]['nHe']
                    next_xval=alldata[next_intermediate][enstring][0]+alldata[next_intermediate]['nHe']

                if current_xval < lowest_xval:
                    lowest_xval=current_xval

                #print(current_intermediate,current_xval)
#                plt.plot([istep+0.75,istep+1.25],
                plt.plot([current_xval-0.25,current_xval+0.25],
                        [En_IS+IS_normalization,En_IS+IS_normalization],linestyle=linestyle,
                        color=color,linewidth=4)
                all_ens.append(En_IS+IS_normalization)
                if istep == len(single_mech_steps)-1:
                    continue

                all_ens.append(En_FS+IS_normalization)

                #If the barrer to the FS has not been calculated
                #Draw a straight line to the next intermediate
                arrow_xlen=0.5+(next_xval-current_xval-1)
                straight_linestyle=linestyle
#                print(current_intermediate,barenstring,alldata[current_intermediate].keys())
#                print(barenstring,alldata[current_intermediate][barenstring])
                Barrier_found=0
                if linestyle=='--': straight_linestyle=':' #Hack because -- looks like a solid line
                if barenstring in alldata[current_intermediate]:
                 if next_intermediate in alldata[current_intermediate][barenstring]:# or next_intermediate == 'H':
                  if Eddag is not None:
                    Barrier_found=1
                    x_IS=current_xval+0.25
                    x_FS=next_xval-0.25
                    barx=np.sqrt((En_IS-Eddag)*(En_FS-Eddag)*(x_IS-x_FS)**2)+En_IS*x_FS-En_FS*x_IS+Eddag*(x_IS-x_FS)
                    barx/=En_IS-En_FS
                    if np.isnan(barx):
                        txtout=f'WARNING! Automatic detection of barrier x '
                        txtout+=f'for {current_intermediate} to '
                        txtout+=f'{next_intermediate} at {potential}V '
                        txtout+=f'failed, likely you are '
                        txtout+=f'barrierless. Check with check_barriers_with_line keyword.'
                        print(txtout)
                        barx=current_xval+(next_xval-current_xval)/2+0.08*(En_FS-En_IS)
                    #if add_potential_response_in_x:
                    # x_IS=current_xval+0.25
                    # x_FS=next_xval-0.25
                    # barx=np.sqrt((En_IS-Eddag)*(En_FS-Eddag)*(x_IS-x_FS)**2)+En_IS*x_FS-En_FS*x_IS+Eddag*(x_IS-x_FS)
                    # barx/=En_IS-En_FS
                    parpts=np.array([[current_xval+0.25,En_IS],
                     #   #[current_xval+(next_xval-current_xval)/2+0.04*(En_FS-En_IS),Eddag], # x-value should b/(2a) of the parabola function
                        [barx,Eddag], # x-value should b/(2a) of the parabola function
                        [next_xval-0.25,En_FS]])
                    #else:
                    # barx=np.sqrt((En_IS-Eddag)*(En_FS-Eddag)*(x_IS-x_FS)**2)+En_IS*x_FS-En_FS*x_IS+Eddag*(x_IS-x_FS)
                    # barx/=En_IS-En_FS
                    # if np.isnan(barx):
                     #    barx=current_xval+(next_xval-current_xval)/2+0.08*(En_FS-En_IS)
                     #print(barx)
                     #parpts=np.array([[current_xval+0.25,En_IS],
                     #   [barx,Eddag],
#                   #     [current_xval+(next_xval-current_xval)/2+0.08*(En_FS-En_IS),Eddag], # x-value should b/(2a) of the parabola function
                     #   [next_xval-0.25,En_FS]])

                    from general_tools import quad_fun
                    from scipy.optimize import curve_fit
                    coeff,dummy=curve_fit(quad_fun,parpts[:,0],parpts[:,1])
                    fitpts=np.linspace(current_xval+0.25,next_xval-0.25,30)#,include_endpoints=True)
                    fit=np.poly1d(coeff)(fitpts)+IS_normalization

                    # Emphasize chosen barriers based on emphasize_barrier input
                    emphasize=False
                    linewidth=1
                    if len(emphasize_barriers):
                        if isinstance(emphasize_barriers,dict):
                            for prop in emphasize_barriers:
                                if prop[:2] == 'pH':
                                    if abs(ph - float(prop[2:])) < 0.01:
                                        pairs=emphasize_barriers[prop]
                                else:
                                    print('pH given for emphasizing barriers not understood')
                        elif isinstance(emphasize_barriers,list):
                            pairs=emphasize_barriers

                        for pair in pairs:
                                if current_intermediate in pair and next_intermediate in pair:
                                    linewidth=6

                    # Plot the parabolic barriers
                    # XXX color of acidic needs to be updated
                    if pdonor == 'acid':
                            plt.plot(fitpts,fit,'--b',linewidth=linewidth)
#                            plt.plot([istep+1.35,istep+1.65],[Eddag+IS_normalization,Eddag+IS_normalization],'--b')
                    else:
                            #plt.plot(fitpts,fit,linestyles[imech],color=color)
                            plt.plot(fitpts,fit,linestyle,color=color,linewidth=linewidth)
                            if check_barriers_with_line:
                               plt.plot([next_xval-0.65,next_xval+0.35],[Eddag+IS_normalization,Eddag+IS_normalization],'--k') #+colors[imech])

                    all_ens.append(Eddag+IS_normalization)

                if not Barrier_found:
                    print(f'Barrier for step {current_intermediate} to {next_intermediate} in {pdonor} not found')
                    # Awful hack regarding the color for SelectCO2 deliverable
                    linecolor=colors[imech]
                    if pdonor == 'acid':
                            linecolor='b'
                    plt.arrow(current_xval+0.25,En_IS+IS_normalization,arrow_xlen,
                                En_FS-En_IS,head_width=0.0,length_includes_head=True,
                                linewidth=1.5,color=linecolor,linestyle=straight_linestyle)



    plt.ylim(ylim)
    if not add_potential_response_in_x:
        plt.xlim([-0.25,longest_mech-lowest_xval-0.75])
    else:
        #print(longest_mech-lowest_xval)
        plt.xlim([lowest_xval-0.5,longest_mech-lowest_xval+2.5])
        plt.xlabel(xlabel)


    if xticks is None:
        plt.xticks([])
    else:
        plt.xticks(np.unique(_xticks))
    if proton_donor=='base':
        donor='H$_2$O'
    elif proton_donor == 'acid':
        donor='H$_3$O$^+$'
    elif proton_donor == 'chemical':
        donor='chemical'

    elif isinstance(proton_donor,dict):
        donor=list(proton_donor.items())
#    print(','.join(['-'.join(i) for i in included_steps]))
    if title is None:
        plt.title(#','.join(['-'.join(i) for i in included_steps])+
            ', facet: '+facet+
            ', WF='+'-'.join([str(np.around(i,2)) for i in potentials])+
            ', pH='+'-'.join([str(np.around(i,1)) for i in pH])+
            ', proton donor: %s '%donor,fontsize=15)
        #'-'.join(included_steps))+
    else:
        plt.title(title)

    plt.ylabel(f'$\Delta${energy_type}$^\phi$ / eV')

    if annotate_reaction_conditions:
        sheout=','.join([str(np.around(i-V0_SHE,2)) for i in potentials])+'V$_{\mathrm{SHE}}$'
        phout = ','.join([str(np.around(i,1)) for i in pH])
        plt.annotate(f"Cu({facet})\n{sheout}\npH={phout}",(-0.2,ylim[0]),ha='left',va='bottom',fontsize=23).draggable()

    plt.tight_layout()
    if view:
        if return_plt:
            return plt
        else:
            plt.show()
            return
    if potentials:
        #print(included_steps)
        if isinstance(included_steps[0],list):
            stepsout=['-'.join(steps) for steps in included_steps]
        else:
            stepsout=included_steps
        #print(stepsout)
        plt.savefig(outdir+'/FED_w_barrier_'+
                '_'.join(stepsout)+
                '_pot_'+'-'.join([str(np.around(i,2)) for i in potentials])+
                '_pH_'+'-'.join([str(np.around(i,2)) for i in pH])+
                '.'+outformat,transparent=True)
    else:
        plt.savefig(outdir+'/FED_w_barrier_'+'-'.join(included_steps)+'.'+outformat,transparent=True)
    plt.close()

#main()
def read_beta_from_catmap_input(line):
            beta=None
            if 'beta' in line:
                beta_string=line.split('beta=(')[1]
                beta_string=beta_string.split(')')[0]
                beta=[float(i.replace(',','')) for i in beta_string.split()]
            return beta

