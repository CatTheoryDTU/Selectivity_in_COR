from ase.io import read,write
from scipy.optimize import curve_fit,OptimizeWarning

import matplotlib
import pickle
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import os,sys
import numpy as np
from .random_utilities import get_reference_energies,lin_fun,quad_fun,get_reference_vibrational_contribution
from ase import units
import warnings

#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (7,5)
markersize=10

from matplotlib.ticker import FormatStrFormatter

intcol={'OHdes':'y','Oprot':'b','Cprot':'brown','desorption':'g','else':'k','unknown':'r'}
def plot_barrier_vs_adsorption_energy(alldata,facet,outfile,specific_ads=[]):
    return None
    colors=cm.gist_ncar(np.linspace(0,1,len(alldata.keys())+1))
    for i,ads in enumerate(alldata.keys()):
        if specific_ads:
            if 'clean' not in specific_ads:
                specific_ads.append('clean')
            if ads not in specific_ads:
                continue



def plot_E_vs_q(alldata,facet,outfile,transform_to_capacity=True,plot_relative_to_clean=False,specific_ads=[],ylim=None):

    colors=cm.gist_ncar(np.linspace(0,1,len(alldata.keys())+1))
    Surf_area = np.product(np.diag(alldata['clean']['cell'][:2,:2]))
    for i,ads in enumerate(alldata.keys()):
        if ads == 'clean':
            continue
        if 'ne_%s'%facet not in alldata[ads].keys():
            continue
        if specific_ads:
            if 'clean' not in specific_ads:
                specific_ads.append('clean')
            if ads not in specific_ads:
                continue
        Cap_data=[]
        for pot in alldata[ads]['ne_%s'%facet].keys():
            if transform_to_capacity:
                Cap_data.append([alldata[ads]['ne_%s'%facet][pot]*-1.6022*1e3/Surf_area,alldata[ads]['E_C_%s'%facet][alldata[ads]['ne_%s'%facet][pot]]])
        Cap_data=np.array(Cap_data)
        if len(Cap_data) < 2: continue
        coeff,d = curve_fit(lin_fun,Cap_data[:,0],Cap_data[:,1])
        fit=np.array([[min(Cap_data[:,0]),min(Cap_data[:,0])*coeff[0]+coeff[1]],
            [max(Cap_data[:,0]),max(Cap_data[:,0])*coeff[0]+coeff[1]]])
        plt.plot(fit[:,0],fit[:,1],'-',color=colors[i])
        if transform_to_capacity:
            plt.plot(Cap_data[:,0],Cap_data[:,1],'+',color=colors[i],label=ads+r', %1.2f $eV/(\mu C/cm^2)$'%(coeff[0]))
        elif plot_relative_to_clean:
            plt.plot(Cap_data[:,0],Cap_data[:,1],'+',color=colors[i],label=ads+r',$ \partial ne_{rel} / \partial\Phi$=%1.2f'%coeff[0])
        else:
            plt.plot(Cap_data[:,0],Cap_data[:,1],'+',color=colors[i],label=ads)

    if ylim:
        plt.ylim(ylim)
    plt.xlabel(r'Surface charge [$\mu$C/cm$^2$]')
    if transform_to_capacity:
        plt.ylabel(r'Formation energy [eV]')
    elif plot_relative_to_clean:
        plt.ylabel(r'Excess electrons relative to clean [e]')
    else:
        plt.ylabel(r'Excess electrons in cell [e]')
    plt.legend(loc='lower left',bbox_to_anchor=(1, 0),fontsize=6)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def plot_cavity(alldata,facet,outfile,specific_ads=[],highlight_ads=[],linewidth=0.5):
    colors=cm.gist_ncar(np.linspace(0,1,len(alldata.keys())+1))

    highlight_ads.append('clean')
    for i,ads in enumerate(alldata.keys()):

        if 'cavity_%s'%facet not in alldata[ads].keys():
            continue

        if specific_ads:
            if 'clean' not in specific_ads:
                specific_ads.append('clean')
            if ads not in specific_ads:
                continue
        if ads not in highlight_ads:
            plt.plot(alldata[ads]['cavity_%s'%facet][:,0],alldata[ads]['cavity_%s'%facet][:,1],'-',linewidth=linewidth,label=ads,color=colors[i])

    for i, ads in enumerate(highlight_ads):
        if 'cavity_%s'%facet not in alldata[ads].keys():
            print('cavity data of %s not found!'%ads)
            continue
        plt.plot(alldata[ads]['cavity_%s'%facet][:,0],alldata[ads]['cavity_%s'%facet][:,1],'-k',linewidth=linewidth*2,label=ads)

    plt.legend(loc='lower left',bbox_to_anchor=(1, 0),fontsize=6)
    plt.xlim([0,alldata['clean']['cell'][2,2]])
    plt.ylabel('Cavity shape function')
    plt.xlabel(r'z-direction [$\AA{}$]')
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def plot_capacity(alldata,facet,outfile,transform_to_capacity=True,plot_relative_to_clean=False,specific_ads=[],ylim=None,fontsize=4):

    plt.close()
    colors=cm.gist_ncar(np.linspace(0,1,len(alldata.keys())+1))
    Surf_area = np.product(np.diag(alldata['clean']['cell'][:2,:2]))
    for i,ads in enumerate(alldata.keys()):

        if 'ne_%s'%facet not in alldata[ads].keys():
            continue
        if specific_ads:
            if 'clean' not in specific_ads:
                specific_ads.append('clean')
            if ads not in specific_ads:
                continue
        Cap_data=[]
        for pot in alldata[ads]['ne_%s'%facet].keys():
            if transform_to_capacity:
                Cap_data.append([pot, alldata[ads]['ne_%s'%facet][pot]*-1.6022*1e3/Surf_area])
            elif plot_relative_to_clean:
                Cap_data.append([pot, alldata[ads]['ne_%s'%facet][pot]-alldata['clean']['ne_%s'%facet][pot]])
            else:
                Cap_data.append([pot, alldata[ads]['ne_%s'%facet][pot]])
        Cap_data=np.array(Cap_data)
        coeff,d = curve_fit(lin_fun,Cap_data[:,0],Cap_data[:,1])
        fit=np.array([[min(Cap_data[:,0]),min(Cap_data[:,0])*coeff[0]+coeff[1]],
            [max(Cap_data[:,0]),max(Cap_data[:,0])*coeff[0]+coeff[1]]])
        plt.plot(fit[:,0],fit[:,1],'-',color=colors[i])
        if transform_to_capacity:
            plt.plot(Cap_data[:,0],Cap_data[:,1],'+',color=colors[i],label=ads+r', %1.2f $\mu F/cm^2$'%coeff[0])
        elif plot_relative_to_clean:
            plt.plot(Cap_data[:,0],Cap_data[:,1],'+',color=colors[i],label=ads+r',$ \partial ne_{rel} / \partial\Phi$=%1.2f'%coeff[0])
        else:
            plt.plot(Cap_data[:,0],Cap_data[:,1],'+',color=colors[i],label=ads)
    plt.xlim([2.4,3.6])
    if ylim:
        plt.ylim(ylim)
    plt.xlabel('Work function [eV]')
    if transform_to_capacity:
        plt.ylabel(r'Surface charge density [$\mu$C/cm$^2$]')
    elif plot_relative_to_clean:
        plt.ylabel(r'Excess electrons relative to clean [e]')
    else:
        plt.ylabel(r'Excess electrons in cell [e]')
    plt.legend(loc='lower left',bbox_to_anchor=(1, 0),fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def get_E_vs_pot(alldata,ads,potentials,enkey,deeperkey=None):
        E_v_pot=[]
        for potential in np.around(potentials,2):
            if deeperkey is not None:
              if isinstance(deeperkey,str):
                if str(potential) in alldata[ads][enkey][deeperkey].keys():
                    E_v_pot.append([potential,alldata[ads][enkey][deeperkey][str(potential)]])
              #if the dict goes deeper (e.g. for barriers)
              elif isinstance(deeperkey,list):
                    fulldepth = len(deeperkey)
                    en=alldata[ads][enkey]
                    for depth in range(fulldepth):
                        en=en[deeperkey[depth]]
                    if str(potential) in en.keys():
                        E_v_pot.append([potential,en[str(potential)]])
                    #This is for accounting for directory names "pot_2.40" instead of "pot_2.4"
                    elif '%1.2f'%potential in en.keys():
                        E_v_pot.append([potential,en['%1.2f'%potential]])

            else:
                if potential in alldata[ads][enkey].keys():
                    #print(alldata[ads][enkey][potential],enkey)
                    #das
                    E_v_pot.append([potential,alldata[ads][enkey][potential]])
        return np.array(E_v_pot)

def fit_potential_response(alldata,facet,potentials,include_barriers=True,plot=True,plotoutname='E_',specific_ads=None,quad_fit=False,not_enough_potentials=None):
    #TURN THIS THE FOLLOWING OFF FOR DEBUGGING!
    warnings.simplefilter('ignore',category=OptimizeWarning)
    if not_enough_potentials is None:
        not_enough_potentials=[]
    for iads,ads in enumerate(alldata.keys()):
     #If only a single adsorbate should be fitted
     if specific_ads:
         if isinstance(specific_ads,str):
             specific_ads=[specific_ads]
         if ads not in specific_ads: continue
     if ads[-2:] == '_g': continue
     if 'E_%s'%facet in alldata[ads].keys():
        E_v_pot = get_E_vs_pot(alldata,ads,potentials,'E_%s'%facet)
        if quad_fit:
            if len(E_v_pot) > 2:
                coeff,dummy = curve_fit(quad_fun,E_v_pot[:,0],E_v_pot[:,1])
                alldata[ads]['E_vs_pot_%s'%facet]=coeff
                if len(E_v_pot) == 3:
                    not_enough_potentials.append(ads)
        else:
            if len(E_v_pot) > 1:
                coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
                alldata[ads]['E_vs_pot_%s'%facet]=coeff
            if len(E_v_pot) == 2:
                not_enough_potentials.append(ads)

     if 'E_C_%s'%facet in alldata[ads].keys():
        E_v_pot = get_E_vs_pot(alldata,ads,potentials,'E_%s'%facet)
        if len(E_v_pot) > 1:
            coeff,onset = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
            alldata[ads]['E_C_vs_pot_%s'%facet]=coeff
        if len(E_v_pot) == 2:
                not_enough_potentials.append(ads)
    #if len(not_enough_potentials):
    #    print('More potentials should be calcualted for: ', not_enough_potentials)

    if include_barriers:
        alldata = fit_barriers_vs_pot(alldata,facet,potentials,plot=plot)

    if plot:
        plot_E_vs_pot(alldata,'E',facet,potentials,quad_fit=quad_fit)
        plot_betas(alldata,facet)

    #warnings.simplefilter('default')#,category='OptimizeWarning')
    return alldata,not_enough_potentials

def plot_betas(alldata,facet=None,exclude=[],specific_reactions=[],
        outdir='results/',ylabel='Potential response / (eV/V)',ax=None,return_plt=False,
        add_legend=True,markers=None,xticks=None,plot_meanval_by_type=False,meanval_label=True,
        exclude_reaction_types=[],markersize=8,plot_thermodynamics=True,plot_kinetics=True):

    if facet is None:
        facet='100'

    #GetThermo beta
    if plot_thermodynamics:
        xlabels=[]
        betas=[]
        for ads in alldata.keys():
            if ads in exclude: continue
            if ('E_vs_pot_%s'%facet not in alldata[ads].keys() or
            ads in ['clean']): continue
            ads_short = ads.lstrip('md-').lstrip('bdo-')
            betas.append(alldata[ads]['E_vs_pot_%s'%facet][0])
            xlabels.append(ads_short)


        fig_thermo,ax_thermo = plt.subplots(figsize=(len(xlabels)/2.,len(xlabels)/4.))

        sortbeta = np.argsort(np.array(betas))
        xlabels = [xlabels[i] for i in sortbeta]
        betas = [betas[i] for i in sortbeta]
        for ibeta,beta in enumerate(betas):
            ax_thermo.plot(ibeta,beta,'o',
                 markeredgecolor='k',color='k',markersize=8)

        ax_thermo.set_xlim(-0.5,len(xlabels)-0.5)
        ax_thermo.set_xticklabels(ax_thermo.get_xticks(),rotation=90)
        ax_thermo.set(xticks=np.arange(len(xlabels)))
        ax_thermo.set(xticklabels=xlabels)
        ax_thermo.axhline(y=0,linestyle='--',color='k')

        ax_thermo.set_ylabel(ylabel,fontsize=16)
        #plt.legend()
        if not return_plt:
            plt.tight_layout()
            if facet is None:
                fig_thermo.savefig(f'{outdir}/Beta_thermo.pdf')
            else:
                fig_thermo.savefig(f'{outdir}/Beta_thermo_%s.pdf'%facet)
            plt.close()
        elif not plot_kinetics:
            return ax_thermo

    if ax is None:
        fig_kin,ax_kin = plt.subplots()
    else:
        ax_kin=ax

    ax_kin.plot(np.nan,np.nan,'*',color=get_intcolors('CO','CHO'),label='C protonation')
    ax_kin.plot(np.nan,np.nan,'*',color=get_intcolors('COH','C'),label='OH desorption')
    ax_kin.plot(np.nan,np.nan,'*',color=get_intcolors('CO','COH'),label='O protonation')
    ax_kin.plot(np.nan,np.nan,'*',color=get_intcolors('HCCO','H2CCO'),label='Anionic desorption')

    xlabels=[]
    betas=[]
    colors,markersout=[],[]

    beta_by_type={}
    for ads in alldata.keys():
        if 'E_ddag_vs_pot_%s'%facet not in alldata[ads].keys(): continue
        ads_short = ads.lstrip('md-').lstrip('bdo-')
        for toads in alldata[ads]['E_ddag_vs_pot_%s'%facet].keys():
            if get_intcolors(ads,toads,return_name=True) in exclude: continue
            if ads_short == 'H2CCH2O' and toads == 'O': continue


            toads_short = toads.lstrip('md-').lstrip('bdo-')
            #print(markers[ads],ads,toads)
            if len(exclude_reaction_types):
               if get_intcolors(ads,toads) in exclude_reaction_types: continue

            if len(specific_reactions):
                neglect=True
                for spec_rea in specific_reactions:
                    if (ads_short == spec_rea.split('-')[0] and
                        toads_short == spec_rea.split('-')[1]):
                        neglect=False
                if neglect:
                    continue

            if len(exclude):
                neglect=False
                for exc_rea in exclude:
                    if '-' not in exc_rea: continue
                    if (ads_short == exc_rea.split('-')[0] and
                        toads_short == exc_rea.split('-')[1]):
                        neglect=True
                if neglect:    continue
            ##TAKE THIS OUT AFTER!!
            #if get_intcolors(ads,toads) in ['r','k']:
            #    continue

            betas.append(alldata[ads]['E_ddag_vs_pot_%s'%facet][toads]['base'][0])

            if get_intcolors(ads,toads,return_type=True) not in beta_by_type:
                beta_by_type[get_intcolors(ads,toads,return_type=True)]=[]
            beta_by_type[get_intcolors(ads,toads,return_type=True)].append(alldata[ads]['E_ddag_vs_pot_%s'%facet][toads]['base'][0])

            colors.append(get_intcolors(ads,toads))
            if markers is None: markersout.append('o')
            else: markersout.append(markers[ads][toads_short])
            #print(ads,toads,markers[ads])
            xlabels.append(get_intcolors(ads,toads,return_name=True))#'%s-%s'%(ads_short,toads_short))

    meanvals_by_type={}
    for typ in beta_by_type:
        meanvals_by_type[typ] = np.mean(beta_by_type[typ])

    if plot_meanval_by_type:
        typout={'Oprot':['Mean O-H',2],'Cprot': ['Mean C-H',23],'OHdes': ['Mean -OH',18],'desorption': ['',15]}
        for typ in meanvals_by_type:
            if len(beta_by_type[typ]) > 4:
                ax_kin.axhline(y=meanvals_by_type[typ],color=intcol[typ],linestyle=':',linewidth=2)
                if meanval_label:
                    ax_kin.annotate(typout[typ][0],(typout[typ][1],meanvals_by_type[typ]),
                            color=intcol[typ],bbox={'edgecolor':'w','facecolor':'w','pad':0.0},
                            ha='center',va='center',fontsize=15).draggable()

    sortbeta = np.argsort(np.array(betas))
    xlabels = [xlabels[i] for i in sortbeta]
    colors = [colors[i] for i in sortbeta]
    markersout = [markersout[i] for i in sortbeta]
    betas = [betas[i] for i in sortbeta]
    for ibeta,beta in enumerate(betas):

        ax_kin.plot(ibeta,beta,markersout[ibeta],
             markeredgecolor='k',color=colors[ibeta],markersize=markersize)

    if xticks is None:
        ax_kin.set_xticklabels(ax_thermo.get_xticks(),rotation=90)
        ax_kin.set(xticks=np.arange(len(xlabels)))
        ax_kin.set(xticklabels=xlabels)
    else:
        ax_kin.set(xticks=xticks)

#    ax_kin.axhline(y=0.5,linestyle='--',color='k')
    ax_kin.set_ylabel(ylabel,fontsize=18)
    if add_legend: ax_kin.legend(loc='lower right')
    ax_kin.set_ylim([0,1])
    if not return_plt:
        plt.tight_layout()
        if facet is None:
            plt.savefig(f'{outdir}/Beta_barriers.pdf')
        else:
            plt.savefig(f'{outdir}/Beta_barriers_%s.pdf'%facet)
    elif not plot_thermodynamics:
        return ax_kin
    else:
        return ax_thermo,ax_kin
    plt.close()

def fit_barriers_vs_pot(alldata,facet,potentials,plot=False):
    for iads,ads in enumerate(alldata.keys()):
     if 'E_ddag_%s'%facet not in alldata[ads].keys():
         continue
     for outname in ['E_ddag_%s'%facet,'E_ddag_IS_%s'%facet,'E_ddag_FS_%s'%facet]:
         for toads in alldata[ads][outname].keys():

          for pH in ['base','acid','chemical']:
            for quad_fit in [False,True]:

              if pH not in alldata[ads][outname][toads].keys():
                  continue
              en = alldata[ads][outname][toads][pH]
              if len(en.keys()) > 2:
                 E_v_pot = get_E_vs_pot(alldata,ads,potentials,outname,deeperkey=[toads,pH])
                 #print(ads,pH,en,E_v_pot)
                 if quad_fit:
                     if len(E_v_pot) < 3:
                         print("Quadratic fit for %s to %s failed due to lack of  potentials"%(ads,toads))
                         continue
                     coeff,dummy = curve_fit(quad_fun,E_v_pot[:,0],E_v_pot[:,1])
                 else:
                     coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
              elif quad_fit: continue
              elif len(en.keys()) > 1:
                 E_v_pot = get_E_vs_pot(alldata,ads,potentials,outname,deeperkey=[toads,pH])
                 if len(E_v_pot) > 1:
                     coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
                 else:
                     continue
              else:
                  continue

              fitname=outname.replace(facet,'')+'vs_pot'
              if quad_fit:
                  fitname+='_quad'#E_ddag_vs_pot_quad'
              fitname+='_%s'%facet
              if fitname not in alldata[ads].keys():
                     alldata[ads][fitname]={}
              if toads not in alldata[ads][fitname].keys():
                     alldata[ads][fitname][toads]={}
              alldata[ads][fitname][toads][pH]=coeff
              if 'E_vs_pot_%s'%facet not in alldata[ads]:
                  print('Thermodynamics of IS for %s to %s on facet %s are missing'
                          %(ads,toads,facet))
                  continue

              E_v_pot2=np.array(E_v_pot)
              E_v_pot2[:,1]=1
              E_v_pot[:,1]-=E_v_pot2@np.array(alldata[ads]['E_vs_pot_%s'%facet])

              relname=outname.replace(facet,'')+'rel_%s'%facet
              if relname not in alldata[ads].keys():
                     alldata[ads][relname]={}
              if toads not in alldata[ads][relname].keys():
                     alldata[ads][relname][toads]={}
              alldata[ads][relname][toads][pH]={}
              for pot_E in E_v_pot:
                  alldata[ads][relname][toads][pH]['%s'%pot_E[0]]=pot_E[1]

              fitname=outname.replace(facet,'')+'rel_vs_pot'
              if quad_fit:
                  fitname+='_quad'#E_ddag_vs_pot_quad'
              fitname+='_%s'%facet
              if quad_fit:
                coeff,dummy = curve_fit(quad_fun,E_v_pot[:,0],E_v_pot[:,1])
              else:
                  coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
              if fitname not in alldata[ads].keys():
                     alldata[ads][fitname]={}
              if toads not in alldata[ads][fitname].keys():
                     alldata[ads][fitname][toads]={}
              alldata[ads][fitname][toads][pH]=coeff
    if plot:
        quad_fit=False
        plot_E_vs_pot(alldata,'E_ddag',facet,potentials,True,quad_fit)
        plot_E_vs_pot(alldata,'E_ddag_rel',facet,potentials,True,quad_fit,ylabel='E$_a$ [eV]')
        plot_Eddag_from_specific_IS(alldata,'E_ddag',facet,potentials,True,quad_fit)

    return alldata

    #Plot BEP
def plot_BEP_relation(alldata,facet,output_potential=3.63,pH=14,exclude=[],outdir='results/',
        xlabel=None,ylabel=None,proton_donors=['base','acid','chemical'],return_plt=False,figax=None,plot_legend=True,
        excluded_types=['k','r','g'],annotate_fit_function=False,markersize=8):


    if isinstance(proton_donors,str): proton_donors=[proton_donors]

    if figax is None:
        fig,ax=plt.subplots()
    else:
        fig,ax=figax

    markers=['d', '^', 'h', 'H','>','<','o','8','X','D','3']
    counters=[0,0,0,0,0]

    #For legend of plot
    ax.plot(np.nan,np.nan,'o', color='brown',label='C protonation')
    ax.plot(np.nan,np.nan,'oy',label='OH desorption')
    ax.plot(np.nan,np.nan,'ob',label='O protonation')
    #ax.plot(np.nan,np.nan,'ok',label='Chemical')
    ax.plot(np.nan,np.nan,'ow',label='-'*20+'\n')

    RHE=output_potential+0.059*pH-4.40

    specific_markers={}
    for  protdon in proton_donors:
        dE_dagg_vs_IS,dE_dagg_rel_vs_IS=[],[]
        dE_dagg_vs_FS,dE_dagg_rel_vs_FS=[],[]
        bars_by_type={'Cprot':[],'OHdes':[],'Oprot':[]}
        dE_dagg_vs_Erxn=[]
        names,markerout=[],[]
        for iads,ads in enumerate(alldata.keys()):
         if ads in ['clean']: continue
#         print(alldata[ads].keys(),'E_ddag_%s'%facet)
#         print(alldata[ads]['E_ddag_%s'%facet])
         if 'E_ddag_%s'%facet not in alldata[ads].keys():
             continue

         for toads_in in alldata[ads]['E_ddag_%s'%facet].keys():
          if protdon not in alldata[ads]['E_ddag_%s'%facet][toads_in]: continue

    #      if ads=='HOCCOH' and toads=='CCOH':
    #          print(alldata[ads]['E_ddag_%s'%facet][toads_in][protdon])
#              das
          if len(alldata[ads]['E_ddag_%s'%facet][toads_in][protdon].keys()) <  2:
              print('Not enough potentials have been calculated for %s to %s in %s'%(ads,toads_in,protdon))
              continue

          found_toads_name=None
          for nameads in ['md-','bdo-','']:
              if  nameads+toads_in in alldata.keys():
                  found_toads_name=nameads+toads_in
                  break
          #if ads == 'H2CCOH': print(ads,found_toads_name,toads_in)
          if found_toads_name is None:
              print('The product %s has not been found in the intermediates dict'%toads)
              print(alldata.keys())
              continue

          toads=found_toads_name
          ads_short = ads.lstrip('md-').lstrip('bdo-')
          toads_short = toads.lstrip('md-').lstrip('bdo-')
     #     if ads=='OCCOH' and toads=='CCO': ddd
          #if toads in ['CCH','HCCH']: continue
          if toads in ['HCCH']: continue
          if ads in ['CC']: continue

          if len(exclude):
                neglect=False
                for exc_rea in exclude:
                    if '-' not in exc_rea: continue
                    if (ads_short == exc_rea.split('-')[0] and
                        toads_short == exc_rea.split('-')[1]):
                        neglect=True
                if neglect:    continue

          #Exclude chemical barriers from BEP
          for icol,color in enumerate(['b','brown','k','y','g']):
            if get_intcolors(ads,toads) == color:
                if ads not in specific_markers: specific_markers[ads]={}
                #print('a',ads,toads)
                specific_markers[ads][toads_short]=markers[counters[icol]%len(markers)]
                counters[icol]+=1

          if get_intcolors(ads,toads) in excluded_types: continue
          SHE_pot=output_potential
          ads_E=alldata[ads]['E_vs_pot_%s'%facet].copy()
          ads_E[1]=alldata[ads]['E_vs_pot_%s'%facet][0]*SHE_pot+alldata[ads]['E_vs_pot_%s'%facet][1]

          toads_E=alldata[found_toads_name]['E_vs_pot_%s'%facet].copy()
          toads_E[1]=alldata[found_toads_name]['E_vs_pot_%s'%facet][0]*SHE_pot+alldata[toads]['E_vs_pot_%s'%facet][1]

          dE_rxn=toads_E[1]-ads_E[1]
          #Add CHE contribution
          if toads not in ['OCCO']:
            dE_rxn+=RHE
          if ads == 'COCO' and toads == 'OCCOH': continue
          dE_ddag_rel=alldata[ads]['E_ddag_rel_vs_pot_%s'%facet][toads_in][protdon][0]*SHE_pot+\
                      alldata[ads]['E_ddag_rel_vs_pot_%s'%facet][toads_in][protdon][1]
          dE_ddag=alldata[ads]['E_ddag_vs_pot_%s'%facet][toads_in][protdon][0]*SHE_pot+\
                  alldata[ads]['E_ddag_vs_pot_%s'%facet][toads_in][protdon][1]

#          print(specific_markers)
          if toads not in ['O']:
              names.append('%s_%s'%(ads,toads))
              print(ads,toads)
              markerout.append(specific_markers[ads][toads_short])
              dE_dagg_rel_vs_IS.append([ads_E[1],dE_ddag_rel])
              dE_dagg_vs_IS.append([ads_E[1],dE_ddag])
              dE_dagg_rel_vs_FS.append([toads_E[1],dE_ddag_rel])
              dE_dagg_vs_FS.append([toads_E[1],dE_ddag])
              dE_dagg_vs_Erxn.append([dE_rxn,dE_ddag_rel])
              if get_intcolors(ads,toads) == 'b':
                  bars_by_type['Cprot'].append([dE_rxn,dE_ddag_rel])
              elif get_intcolors(ads,toads) == 'y':
                  bars_by_type['OHdes'].append([dE_rxn,dE_ddag_rel])
              elif get_intcolors(ads,toads) == 'brown':
                  bars_by_type['Oprot'].append([dE_rxn,dE_ddag_rel])
        BEP_fit={}

        #BEP_fit={'Cprot':[0,0],'OHdes':[0,0],'Oprot':[0,0]}
        for bartype in bars_by_type:
            bars_by_type[bartype] = np.array(bars_by_type[bartype])
            if len(bars_by_type[bartype]):
                BEP_fit[bartype],dummy = curve_fit(lin_fun,bars_by_type[bartype][:,0],bars_by_type[bartype][:,1])

        #if return_plt:
        #    return plt

        _finalize_BEP_plots(dE_dagg_vs_Erxn,names,markerout,
                outfile=f'{outdir}/BEP_%1.1fV_%s.pdf'%(RHE+0.01,protdon),RHE=RHE,pH=pH,BEP_fit=BEP_fit,xlabel=xlabel,ylabel=ylabel,return_plt=return_plt,figax=(fig,ax),markersize=markersize,
                plot_legend=plot_legend,annotate_fit_function=annotate_fit_function)

        if return_plt: return ax,specific_markers
        if protdon == 'base':
             _finalize_BEP_plots(dE_dagg_rel_vs_IS,names,markerout,xlabel=r'$\Delta$ E$^{IS}_{%1.1fV_{RHE},pH=%i}$'%(RHE,pH),
                outfile=f'{outdir}/E_ddag_rel_vs_E_IS.pdf',RHE=RHE,pH=pH)
             _finalize_BEP_plots(dE_dagg_rel_vs_FS,names,markerout,xlabel=r'$\Delta$ E$^{FS}_{%1.1fV_{RHE},pH=%i}$'%(RHE,pH),
                outfile=f'{outdir}/E_ddag_rel_vs_E_FS.pdf',RHE=RHE,pH=pH)
             _finalize_BEP_plots(dE_dagg_vs_IS,names,markerout,xlabel=r'$\Delta$ E$^{IS}_{%1.1fV_{RHE},pH=%i}$'%(RHE,pH),
                outfile=f'{outdir}/E_ddag_vs_E_IS.pdf',RHE=RHE,pH=pH)
             _finalize_BEP_plots(dE_dagg_vs_FS,names,markerout,xlabel=r'$\Delta$ E$^{FS}_{%1.1fV_{RHE},pH=%i}$'%(RHE,pH),
                outfile=f'{outdir}/E_ddag_vs_E_FS.pdf',RHE=RHE,pH=pH)

def _finalize_BEP_plots(data,names,marker,xlabel=None,ylabel=None,
        outfile='results/BEP.pdf',RHE=0,pH=14,BEP_fit=None,markersize=8,
        return_plt=False,figax=None,plot_legend=True,annotate_fit_function=True):

    if figax is None:
        fig,ax=plt.subplots()
    else:
        fig,ax=figax

    if xlabel is None:
        xlabel=r'$\Delta$ E$^{rxn}_{%1.1fV_{RHE},pH=%i}$'%(RHE+0.01,pH)
    if ylabel is None:
        ylabel=r'$\Delta$ E$^{\dagger}_{%1.1fV_{RHE},pH=%i}$'%(RHE+0.01,pH)

    if BEP_fit:
        colors={'Cprot':'b','Oprot':'brown','OHdes':'y'}
        minx,maxx=min(np.array(data)[:,0]),max(np.array(data)[:,0])
        fitpts=np.linspace(minx,maxx,2)
        for bartype in BEP_fit.keys():
            yvals=BEP_fit[bartype][0]*fitpts+BEP_fit[bartype][1]
            lineangle=np.arctan((yvals[-1]-yvals[0])/(fitpts[-1]-fitpts[0]))*360./(2*np.pi)*0.85
            ax.plot(fitpts,yvals,'--',color=colors[bartype])
            if annotate_fit_function:
                ax.annotate(f'{BEP_fit[bartype][0]:1.2f}$\Delta$E+{BEP_fit[bartype][1]:1.2f}eV' ,
                    (fitpts[-1]-0.1,yvals[-1]),
                    ha='right',va='top',color=colors[bartype],
                    rotation=lineangle,fontsize=16,
                    bbox={'facecolor':'w','edgecolor':'none'}


                    ).draggable()
      #      print(bartype,lineangle)

    for idat,dat in enumerate(data):
        ads,toads=names[idat].split('_')#[0],names[idat].split('-')[1]
        ads_short = names[idat].split('_')[0].lstrip('md-').lstrip('bdo-')
        toads_short = names[idat].split('_')[1].lstrip('md-').lstrip('bdo-')
        ax.plot(dat[0],dat[1],marker[idat],label=get_intcolors(ads,toads,return_name=True),#'%s-%s'%(ads_short,toads_short),
            color=get_intcolors(ads,toads),markersize=markersize,markeredgecolor='k')

    if plot_legend:
        ax.legend(bbox_to_anchor=(1,1.01),fontsize=8.2)
    ax.set_xlabel(xlabel,fontsize=20)
    ax.set_ylabel(ylabel,fontsize=20)

    if return_plt: return ax
    fig.tight_layout()
    #plt.show()
    fig.savefig(outfile)
#    fig.close()

#    plt.plot()

def plot_intrinsic_barrier_energies(alldata,facet,proton_donors=['base','acid','chemical'],
    ylabel='Intrinsic reaction barrier / eV',outdir='results',return_plt=False,figax=None,pH=13,zerovolt=4.40,
    markers=None,xticks=None):

    if figax is None:
        fig,ax=plt.subplots()
    else:
        fig,ax=figax

    for prdon in proton_donors:
     intcounter=0
     intbars={'Cprot':{},'OHdes':{},'Oprot':{}}
     intrbars,colors=[],[]
     xlabels=[]
     markersout=[]
     for iads,ads in enumerate(alldata.keys()):
         if 'E_ddag_%s'%facet not in alldata[ads].keys():
             continue
         if ads in ['clean']: continue
         for toads_in in alldata[ads]['E_ddag_%s'%facet].keys():
          if prdon not in alldata[ads]['E_ddag_%s'%facet][toads_in]:  continue
          if len(alldata[ads]['E_ddag_%s'%facet][toads_in][prdon].keys()) <  2:  continue

          #Add intrinsic barrier
          #Find equilibrium potential
          foundname=None
          for adname in ['md-','bdo-','']:
              if adname+toads_in in alldata.keys():
                  foundname=adname+toads_in
                  break
          if foundname is None:
                  print('Product %s has not been found in intermedates'%toads)
                  continue

          toads=foundname
          SHE_pot=zerovolt-0.059*pH
          ads_E=alldata[ads]['E_vs_pot_%s'%facet].copy()
          ads_E[1]=alldata[ads]['E_vs_pot_%s'%facet][0]*SHE_pot+alldata[ads]['E_vs_pot_%s'%facet][1]
          if ads in ['clean']: ads_E=[0,0]

          toads_E=alldata[toads]['E_vs_pot_%s'%facet].copy()
          toads_E[1]=alldata[toads]['E_vs_pot_%s'%facet][0]*SHE_pot+alldata[toads]['E_vs_pot_%s'%facet][1]

          #Countering CHE for COCO to OCCO:
          if ads not in ['COCO'] and toads not in ['OCCO']:
              toads_E[0]+=1

          E_diff = toads_E-ads_E
          Eq_pot = -E_diff[1]/E_diff[0] + SHE_pot
          if 'Eq_pot_%s'%facet not in alldata[ads].keys():
                 alldata[ads]['Eq_pot_%s'%facet]={}
          alldata[ads]['Eq_pot_%s'%facet][toads]=Eq_pot

          #Calculate intrinsic barrier
          bartype =  get_intcolors(ads,toads)
          if toads_in not in alldata[ads]['E_ddag_rel_vs_pot_%s'%facet]: continue
          if prdon not in alldata[ads]['E_ddag_rel_vs_pot_%s'%facet][toads_in]: continue

          if 'E_ddag_rel_intrinsic_%s'%facet not in alldata[ads].keys():
                 alldata[ads]['E_ddag_rel_intrinsic_%s'%facet]={}
          if toads not in alldata[ads]['E_ddag_rel_intrinsic_%s'%facet].keys():
                 alldata[ads]['E_ddag_rel_intrinsic_%s'%facet][toads_in]={}

          ads_short = ads.lstrip('md-').lstrip('bdo-')
          toads_short = toads.lstrip('md-').lstrip('bdo-')

          alldata[ads]['E_ddag_rel_intrinsic_%s'%facet][toads_in][prdon]=\
                alldata[ads]['E_ddag_rel_vs_pot_%s'%facet][toads_in][prdon][0]*Eq_pot+\
                alldata[ads]['E_ddag_rel_vs_pot_%s'%facet][toads_in][prdon][1]
          for bt in intbars.keys():
              if  prdon not  in  intbars[bt]:
                  intbars[bt][prdon] = []
#          print(ads)
          if bartype == 'brown' and ads not in ['CC','CCH','md-HCCO']: # and toads not in ['H2CCO']:
                  intbars['Cprot'][prdon].append(alldata[ads]['E_ddag_rel_intrinsic_%s'%facet][toads_in][prdon])
          elif bartype == 'y': intbars['OHdes'][prdon].append(alldata[ads]['E_ddag_rel_intrinsic_%s'%facet][toads_in][prdon])
          elif bartype == 'b': intbars['Oprot'][prdon].append(alldata[ads]['E_ddag_rel_intrinsic_%s'%facet][toads_in][prdon])

          #plt.plot(intcounter,alldata[ads]['E_ddag_rel_intrinsic_%s'%facet][toads][pH],
          #          'o',color=get_intcolors(ads,toads),markeredgecolor='k')
          if ads not in ['CC','CCH'] and bartype not in ['k','r','g'] and toads not in ['H2CCO']:
#          if bartype not in ['k']:
              intrbars.append(alldata[ads]['E_ddag_rel_intrinsic_%s'%facet][toads_in][prdon])
              colors.append(get_intcolors(ads,toads))
              markersout.append(markers[ads][toads_short])
              xlabels.append(get_intcolors(ads,toads,return_name=True))#'%s-%s'%(ads_short,toads_short))
              intcounter+=1

     sortbars=np.argsort(np.array(intrbars))
     xlabels = [xlabels[i] for i in sortbars]
     colors = [colors[i] for i in sortbars]
     markersout = [markersout[i] for i in sortbars]
     intrbars = [intrbars[i] for i in sortbars]

     avgcolors=['brown','y','b']
     for ibar,bartype in enumerate(intbars.keys()):
          if prdon not in intbars[bartype]:
              print('No '+prdon+' barriers for  '+bartype+'have been found')
              continue
          mean_intbar=np.mean(intbars[bartype][prdon])
          ax.axhline(y=mean_intbar,linestyle=':',linewidth=2,color=avgcolors[ibar])

          outtxt='Other'
          if bartype=='Cprot': outtxt='Mean C-H'
          elif bartype == 'Oprot': outtxt='Mean O-H'
          elif bartype == 'OHdes': outtxt='Mean -OH'

          ax.annotate(outtxt,(len(xlabels)-1,mean_intbar),ha='right',va='center',color=avgcolors[ibar],bbox={'facecolor':'w','edgecolor':'none'},fontsize=15).draggable()
          print(bartype,mean_intbar)

     for ibar,bar in enumerate(intrbars):
         ax.plot(ibar,bar,'o' if markers is None else markersout[ibar],markersize=8,markeredgecolor='k',color=colors[ibar])
     ax.plot(np.nan,np.nan,'o',color=get_intcolors('COH','C'),label='OH desorption')
     ax.plot(np.nan,np.nan,'o',color=get_intcolors('CO','COH'),label='O protonation')
     ax.plot(np.nan,np.nan,'o',color=get_intcolors('CO','CHO'),label='C protonation')
     #plt.plot(np.nan,np.nan,'o',color=get_intcolors('COCO','OCCO'),label='Other')

     #plt.legend()
     ax.set_ylabel(ylabel,fontsize=18)
     if xticks is None:
         ax.set_xticklabels(ax.get_xticks(),rotation=90)
         ax.set(xticks=np.arange(len(xlabels)))
         ax.set(xticklabels=xlabels)
         #plt.xticks(np.arange(len(xlabels)),xlabels,rotation='vertical')
     else:
         ax.set(xticks=xticks)

     if return_plt:
         return ax
     fig.tight_layout()
     fig.savefig(f'{outdir}/E_ddag_int_%s.pdf'%prdon)
     plt.close()

    return alldata

def get_intcolors(ads,toads,return_name=False,return_type=False):
      ads_short = ads.lstrip('md-').lstrip('bdo-')
      toads_short = toads.lstrip('md-').lstrip('bdo-')
      #XXX Write this in a general manner!!!!
      #intcol={'OHdes':'y','Oprot':'b','Cprot':'brown','desorption':'g','else':'k','unknown':'r'}
      color=None
      if ads == 'CO':
          if  toads == 'CHO': color,name='Cprot','C-HO'
          elif toads == 'COH': color,name='Oprot','CO-H'
      elif ads == 'COH':
          if toads == 'C':  color,name='OHdes','C-OH'
          elif toads == 'OCCOH':  color,name='else','OC-COH'
      elif ads_short == 'CCO':
          if toads_short == 'HCCO':color,name='Cprot','H-CCO'
          if toads_short == 'CCHO':color,name='Cprot','CC-HO'
          elif toads_short == 'CCOH':color,name='Oprot','CCO-H'
      elif ads_short == 'COCO':
          if toads_short == 'OCCO': color,name='else','CO-CO'
      elif ads_short == 'OCCO':
          if toads_short == 'OCCOH':color,name='Oprot','OCCO-H'
      elif ads_short == 'HOCCOH':
          if toads_short=='CCOH': color,name='OHdes','HO-CCOH'
      elif ads_short == 'H2CCH2O':
          if toads_short == 'H3CCH2O': color,name = 'Cprot','H-'+ads_short
          elif toads_short == 'O': color,name = 'else','H2CCH2-O'
      elif ads_short == 'OCCOH':
          if toads_short == 'CCO': color,name = 'OHdes','OCC-OH'
          elif toads_short == 'HOCCOH': color,name = 'Oprot','H-OCCOH'
      elif ads_short == 'CCOH':
          if toads_short == 'CC': color,name = 'OHdes','CC-OH'
          elif toads_short == 'HCCOH': color,name = 'Cprot','H-CCOH'
      elif ads_short == 'clean':
          if toads_short == 'H': color,name = 'else','*-H'
      elif ads_short == 'C':
          if toads_short == 'CH': color,name = 'Cprot','C-H'
      elif ads_short == 'HCCO':
          if toads_short == 'HCCHO': color,name = 'Cprot','HCC-HO'
          elif toads_short == 'HCCOH': color,name = 'Oprot','HCCO-H'
          elif toads_short == 'H2CCO': color,name = 'desorption','H-HCCO'
      elif ads_short == 'CHO':
          if toads_short == 'CHOH': color,name = 'Oprot','CHO-H'
          elif toads_short == 'OCCHO': color,name = 'else','CHO-CO'
      elif ads_short == 'CHOH':
          if toads_short == 'CH': color,name = 'OHdes','CH-OH'
      elif ads_short == 'CC':
          if toads_short == 'CCH': color,name = 'Cprot','CC-H'
      elif ads_short == 'CCH':
          if toads_short == 'HCCH': color,name = 'Cprot','H-CCH'
      elif ads_short == 'OHCCH2':
          if toads_short == 'H2CCH2O': color,name = 'Cprot','H2CCH-HO'
          elif toads_short == 'OCHCH3': color,name = 'desorption','H-CH2CHO'
      elif ads_short == 'H2CCO':
          if toads_short == 'H2CCOH': color,name = 'Oprot','H2CCO-H'
          elif toads_short == 'OHCCH2': color,name = 'Cprot','H2CC-HO'
          elif toads_short == 'H3CCO': color,name = 'Cprot','H-H2CCO'
      elif ads_short == 'H2CCOH':
          if toads_short == 'CCH2': color,name = 'OHdes','H2CC-OH'
      elif ads_short == 'HCCOH':
          if toads_short == 'CCH':  color,name = 'OHdes','HCC-OH'
      else:
          #color,name = 'unknown','CO-CHO'
          color,name='unknown',ads_short+toads_short

      if color is None:
          color,name='unknown',ads_short+toads_short

      if return_name:
          return name
      if return_type:
          return color
      return intcol[color]


def plot_E_vs_pot(alldata,enkey,facet,potentials,deeperkey=None,quad_fit=False,plot_separate=True,ylabel='Formation energy [eV]',specific_adsorbates=None):
    counter=0
    for iads,ads in enumerate(alldata.keys()):
     if enkey+'_%s'%facet in alldata[ads].keys():
         counter+=1

    colors=cm.gist_ncar(np.linspace(0,1,counter+1))
    counter=0
    symbols=['+','s','d']
    for iads,ads in enumerate(alldata.keys()):
     if ads[-2:] == '_g': continue

     if specific_adsorbates is not None:
         if isinstance(specific_adsorbates,str):
             specific_adsorbate=[specific_adsorbates]
         if ads not in specific_adsorbates:
             continue


     if enkey+'_vs_pot_%s'%facet in alldata[ads].keys() and ads not in ['clean']:
        if deeperkey:
          counter2=0
          #print(ads,enkey)
          for toads in alldata[ads][enkey+'_vs_pot_%s'%facet].keys():
           for pH in alldata[ads][enkey+'_vs_pot_%s'%facet][toads].keys():
            coeff=alldata[ads][enkey+'_vs_pot_%s'%facet][toads][pH]
            fit = []
            if quad_fit:
                for pot in np.linspace(potentials[0],potentials[-1],10):
                    fit.append([pot,np.array([pot**2,pot,1])@coeff])
            else:
                for pot in np.linspace(potentials[0],potentials[-1],2):
                    fit.append([pot,np.array([pot,1])@coeff])
            fit=np.array(fit)
            plt.plot(fit[:,0],fit[:,1],color=colors[counter%len(colors)])
            E_v_pot=get_E_vs_pot(alldata,ads,potentials,enkey+'_%s'%facet,[toads,pH])
            if len(E_v_pot):
                plt.plot(E_v_pot[:,0],E_v_pot[:,1],symbols[counter2%len(symbols)],label=ads+'-'+toads+', '+pH+r',dE/d$\Phi$=%1.2f'%coeff[0],color=colors[counter%len(colors)])
                #plt.plot(E_v_pot[:,0],E_v_pot[:,1],symbols[counter%len(symbols)],label=ads+'-'+toads+r',dE/d$\Phi$=%1.2f'%coeff[0],color=colors[counter%len(colors)])
            counter2+=1
        else:
            coeff=alldata[ads][enkey+'_vs_pot_%s'%facet]
            fit = np.array([[potentials[0],potentials[0]*coeff[0]+coeff[1]],
                           [potentials[-1],potentials[-1]*coeff[0]+coeff[1]]])
            plt.plot(fit[:,0],fit[:,1],color=colors[counter%len(colors)])
            E_v_pot=get_E_vs_pot(alldata,ads,potentials,enkey+'_%s'%facet)
            plt.plot(E_v_pot[:,0],E_v_pot[:,1],symbols[counter%len(symbols)],label=ads+r',dE/d$\Phi$=%1.2f'%coeff[0],color=colors[counter%len(colors)])
        counter+=1

    plt.xlabel('Work function [eV]')
    plt.ylabel(ylabel)
    #print(enkey,facet,counter)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4,fontsize=3)
    plt.tight_layout()
    if specific_adsorbates is None:
        plt.savefig('results/%s_vs_pot_%s.pdf'%(enkey,facet))
    else:
        plt.savefig('results/%s_vs_pot_%s.pdf'%(enkey,facet))
    if plot_separate:
        plt.close()


def write_catinput(alldata,facets,potential,outfile,include_barriers=True,
                    catbasein='/Users/geokast/SelectCO2/endstates/tools_for_analysis/catinput_base.txt',
                    catbasefreein='/Users/geokast/SelectCO2/endstates/tools_for_analysis/catinput_base_freeen.txt'):
    out=open(outfile,'w')
    out_G=open('.'.join(outfile.split('.')[:-1])+'_freEn.txt','w')
    basein=open(catbasein,'r').read()
    basefreein=open(catbasefreein,'r').read()
    out.write(basein)
    out_G.write(basefreein)

    if isinstance(facets,str):facets=[facets]
    ##Get the free energy contribution of one water on the clean slab for referencing barriers
    #from ase.thermochemistry import HarmonicThermo
    #from ase.units import invcm

    for facet in facets:
        written_barriers=[]
        for ads in alldata.keys():
            ads_short = ads.lstrip('md-').lstrip('bdo-')
            vib_string=''
            if 'vibs_%s'%facet in alldata[ads].keys():
                vib_string_ads=alldata[ads]['vibs_%s'%facet]
                for vib in alldata[ads]['vibs_%s'%facet]:
                    vib_string+=str(np.around(vib,6))+', '
                vib_string='['+vib_string[:-1]+']'

            if not len(vib_string):
                vib_string = '[]'

            #testprint(ads,facet)
            #try:
            #    testprint(alldata[ads]['E_vs_pot_%s'%facet])
            #except:
            #    testprint(alldata[ads].keys())
            #    sys.exit()
            if 'E_vs_pot_%s'%facet in alldata[ads].keys():

#                testprint(ads,facet)
                coeff = tuple(alldata[ads]['E_vs_pot_%s'%facet])

                out.write('Cu\t%s\t%s\t%f\t%s\tbeta=%s\n'%
                          (facet,ads_short,np.array([potential,1])@coeff,vib_string,coeff))
                if 'free_en_corr_%s'%facet in alldata[ads]:
                    out_G.write('Cu\t%s\t%s\t%f\t%s\tbeta=%s\n'%
                          (facet,ads_short,(np.array([potential,1])@coeff)+alldata[ads]['free_en_corr_%s'%facet]
                              ,'[]',coeff))
                else:
                    print(f'Could not write the G of {ads_short} on {facet} in catmap energy file. missing vibs?')


            if 'E_ddag_vs_pot_%s'%facet in alldata[ads].keys():
                for toads in alldata[ads]['E_ddag_vs_pot_%s'%facet].keys():
                  ads_short = ads.lstrip('md-').lstrip('bdo-')
                  toads_short = toads.lstrip('md-').lstrip('bdo-')

                  for pH in ['base','acid','chemical']:
                    if pH not in alldata[ads]['E_ddag_vs_pot_%s'%facet][toads].keys():
                        continue

                    #if pH == 'acid':
                    #    print('catmap in',ads,toads,
                    #            alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH],
                    #           alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH][0]*4.4+alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH][1])
                    outname,written_barriers,TSnamediff=get_unique_TS_name(ads,toads,pH,written_barriers)
                    coeff=tuple(alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH])
                    out_string='Cu\t%s\t%s\t%f\t%s\tbeta=%s %s_to_%s_%s %s\n'

                    #Add vibrations to transition states
                    ## If no vibrations are found use the IS vibrations for
                    ## and FS vibrations otherwise
                    no_vibs=True
                    if 'vibs_ddag_%s'%facet in alldata[ads]:
                        if toads in alldata[ads]['vibs_ddag_%s'%facet]:
                          if pH in alldata[ads]['vibs_ddag_%s'%facet][toads]:
                            vib_string=alldata[ads]['vibs_ddag_%s'%facet][toads][pH][1:]
                            img_freq=alldata[ads]['vibs_ddag_%s'%facet][toads][pH][0]
                            no_vibs=False

                    if no_vibs:
                            #Add FS vibs as default for acidic barriers
                            if pH in ['acid','chemical']:
                                substitute_endstate = toads
                            #Add IS vibs as default for acidic barriers
                            else:
                                substitute_endstate = ads

                            if 'vibs_%s'%facet in alldata[substitute_endstate]:
                                vib_string=alldata[substitute_endstate]['vibs_%s'%facet]
                            else:
                                vib_string='[]'
#                            else:
#                                vib_string=alldata[ads]['vibs_%s'%facet]
                            img_freq=0

                    if len(coeff) == 3:
                        out.write(out_string%
                              (facet,outname,np.array([potential**2,potential,1])@coeff,
                              vib_string,coeff,ads_short,toads_short,pH,
                              np.around(img_freq,0)))
                    else:
                        out.write(out_string%
                              (facet,outname,np.array([potential,1])@coeff,
                               vib_string,coeff,ads_short,toads_short,pH,
                               np.around(img_freq,0)))
                        #print(ads,out_string)
                    if 'free_en_corr_ddag_%s'%facet in alldata[ads]:
                     #print(ads,toads)
                     if toads in alldata[ads]['free_en_corr_ddag_%s'%facet]:
                        if len(coeff) == 3:
                            out_G.write(out_string%
                                  (facet,outname,np.array([potential**2,potential,1])@coeff+alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH],
                                  '[]',coeff,ads_short,toads_short,pH,
                                  np.around(img_freq,0)))
                        else:
                            out_G.write(out_string%
                                  (facet,outname,np.array([potential,1])@coeff+alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH],
                                  '[]',coeff,ads_short,toads_short,pH,
                                  np.around(img_freq,0)))

def get_unique_TS_name(ads,toads,pH,written_barriers):
    ads_short = ads.lstrip('md-').lstrip('bdo-')
    toads_short = toads.lstrip('md-').lstrip('bdo-')
    ads_short2=ads_short.replace('2','H').replace('3','HH')
    toads_short2=toads_short.replace('2','H').replace('3','HH')

    #Find TS name based on the difference of the names of IS and FS
    for letter in ads_short2:
        toads_short2=toads_short2.replace(letter,'',1)
    outname=ads_short+'-'

    if ads_short == 'clean':
        outname=''

    if len(toads_short2) > 1:
        outname += ads_short
    elif toads_short2 == 'H':
        if pH == 'base':
            outname+='H2O-ele'
        elif pH == 'acid':
            outname+='H-ele'
        elif pH=='chemical':
            outname+='H'
    elif len(toads_short2) == 0:
        if ads_short in ['COCO'] and toads_short in ['OCCO']:
            outname = 'CO-CO'
        else:
            if pH == 'base':
                outname+='ele'
            elif pH == 'acid':
                outname+='H-ele'
            elif pH == 'chemical':
                raise TypeError('Somethings wrong for %s to %s. OH desorption '
                        'can not be chemical')
    else:
        outname+=toads_short

    #Make sure the transition state name is unique
    outname_bakk=outname
    outn=[]
    for i in range(1,100):
        if outname in written_barriers:
            if outname[i] == '-':
                print('Couldnt rename %s barrier from %s to %s'%(pH, ads_short,toads_short))
                outname=outname_bakk
                break
            else:
                outl=list(outname)
                #If the current letter is a number (we don't want to break numbers from letters)
                if outl[i].isnumeric():
                    outn=[j for j in outl[i-1:i+1]]
                    outn+=[letter for j,letter in enumerate(outl)
                          if j not in [i-1,i]]
                #If the following letter is a number
                elif outl[i+1].isnumeric():
                    outn=[j for j in outl[i:i+2]]
                    outn+=[letter for j,letter in enumerate(outl)
                          if j not in [i,i+1]]
                #Else put last letter infront
                else:
                    lastletter=outname.split('-')[0][-1]
                    ilastletter=len(outname.split('-')[0])-1
                    outn=[lastletter]
                    outn+=[letter for j,letter in enumerate(outl)
                          if j != ilastletter]

                outname=''.join(outn)
        else:
            written_barriers.append(outname)
#            if ads == 'md-CCO': print(outname)
            break

    #toads_short2 is the name difference in ads and toads
    return outname,written_barriers,toads_short2

def add_vibrational_free_energy_corrections(alldata,facet,no_water_layer=False,references={'C':'CO'}):
    from ase.thermochemistry import HarmonicThermo
    from ase.units import invcm

    #Get the free energy contribution of one water on the clean slab for referencing barriers
    if not no_water_layer:
        if 'free_en_corr' not in alldata['clean'].keys():
            try:
                vibens=[i*invcm for i in alldata['clean']['vibs_%s'%facet]]
            except:
                vibens=[0]*9#[i*invcm for i in alldata['clean']['vibs_%s'%facet]]
            alldata['clean']['free_en_corr_%s'%facet]=clean_G=HarmonicThermo(vibens).get_helmholtz_energy(298,verbose=False)
        else:
            print('It seems the add_vibrational_free_energy_corrections class is called more than  once')
            return

    #Add vibration correction to adsorbates
    no_freen=[]
    for ads in alldata.keys():
      if ads[-2:] == '_g':
          continue


      if 'E_vs_pot_%s'%facet not in alldata[ads]: continue

      ads_short = ads.lstrip('md-').lstrip('bdo-')
      if ads == 'clean':
          ads_short = ''
      vib_string=''
      vibens=[]
      if 'vibs_%s'%facet in alldata[ads].keys():
          adsvibs=alldata[ads]['vibs_%s'%facet]
          for vib in alldata[ads]['vibs_%s'%facet]:
              vib_string+=str(np.around(vib,4))+', '
              vibens.append(vib*invcm)

          alldata[ads]['free_en_corr_%s'%facet]=HarmonicThermo(vibens).get_helmholtz_energy(298,verbose=False)
          #Subtract gas phase reference free energies
          ads_short_for_vibs=ads_short.replace('2','H').replace('3','HH')
          alldata[ads]['free_en_corr_%s'%facet]-=get_reference_vibrational_contribution(ads_short_for_vibs,references=references)
      else:
          no_freen.append(ads)
          alldata[ads]['free_en_corr_%s'%facet]=0

      alldata[ads]['G_vs_pot_%s'%facet] = alldata[ads]['E_vs_pot_%s'%facet].copy()
      alldata[ads]['G_vs_pot_%s'%facet][1] += alldata[ads]['free_en_corr_%s'%facet]
      if 'E_vs_pot_quad_%s'%facet in alldata[ads]:
          alldata[ads]['G_vs_pot_quad_%s'%facet] = alldata[ads]['E_vs_pot_quad_%s'%facet].copy()
          alldata[ads]['G_vs_pot_quad_%s'%facet][2] += alldata[ads]['free_en_corr_%s'%facet]



      #Transition states
      for ien,enname in enumerate(['_ddag_vs_pot_','_ddag_IS_vs_pot_','_ddag_FS_vs_pot_']):
       Ename,Gname='E'+enname+facet,'G'+enname+facet
       if Ename in alldata[ads].keys():
          for toads in alldata[ads][Ename].keys():
            toads_short = toads.lstrip('md-').lstrip('bdo-')
            for pH in alldata[ads][Ename][toads].keys():
      #        ads_short = ads.lstrip('md-').lstrip('bdo-')
              ## If no vibrations are found use the IS vibrations
              no_vibs=True
              if 'vibs_ddag_%s'%facet in alldata[ads]:
                  if toads in alldata[ads]['vibs_ddag_%s'%facet]:
                    #if pH == 'chemical':
                    #      print(alldata[ads]['vibs_ddag_%s'%facet][toads])
                    if pH in alldata[ads]['vibs_ddag_%s'%facet][toads].keys():
                      vibens=[vib*invcm for vib in alldata[ads]['vibs_ddag_%s'%facet][toads][pH][1:]]
                      no_vibs=False

              #TODO: Add acidic defaulting to FS not IS!
              if no_vibs:
                  vibens=[vib*invcm for vib in adsvibs]
                  if f'vibs_from_IS_{facet}' not in alldata[ads]:
                      alldata[ads][f'vibs_from_IS_{facet}']={}
                  alldata[ads][f'vibs_from_IS_{facet}'][toads]=True
                  no_freen.append(f'{ads}-{toads}{facet}{pH}')
              if not ien: #Only for barrier
                  if 'free_en_corr_ddag_%s'%facet not in alldata[ads]:
                      alldata[ads]['free_en_corr_ddag_%s'%facet]={}
                  if toads not in alldata[ads]['free_en_corr_ddag_%s'%facet]:
                      alldata[ads]['free_en_corr_ddag_%s'%facet][toads]={}
                  eng_corr=HarmonicThermo(vibens).get_helmholtz_energy(298,verbose=False)

                  if pH == 'base':
                      name_for_reference=ads_short.replace('2','H').replace('3','HH')
                  elif pH in ['acid','chemical']:
    #                  print('acid or chemical',ads,toads)
                      name_for_reference=toads_short.replace('2','H').replace('3','HH')

                  outname,dummy,TSnamediff=get_unique_TS_name(ads,toads,pH,[])
                  eng_corr_before_subtraction=eng_corr
                  #The free energy contribution of one water molecule is subtracted, because endstates do not contain
                  #the water, while barriers need to in order to get the free energy
                  if (any(i in outname for i in ['-H2O-','-H-']) or outname == 'H2O-ele') and pH != 'chemical':
                      if not no_vibs: #and (toads not in ['H2CCO'] and not facet == '100'):
                       eng_corr-=clean_G
                      #HCCO to H2CCO is excluded because the desorption not protonation are the barrier
                      if toads in ['H2CCO','OCHCH3'] and facet == '100' and not no_water_layer: #not no_vibs:
                          eng_corr+=clean_G

    #              print(ads,toads,alldata[ads]['free_en_corr_%s'%facet])
                  eng_corr-=get_reference_vibrational_contribution(name_for_reference,references=references)

                  #Add a warning if the free energy correction vary too much
                  if (eng_corr-alldata[ads]['free_en_corr_%s'%facet] < -0.35 and ads != 'clean'):
                      print(f"WARNING! The difference in Free energy correction of the IS and TS of {ads}-{toads} vary by"
                              f" {eng_corr-alldata[ads]['free_en_corr_%s'%facet]:1.2f} (expected 0 > x > 0.3)")
    #              print(ads,toads,eng_corr,alldata[ads]['free_en_corr_%s'%facet],eng_corr-alldata[ads]['free_en_corr_%s'%facet])
                  alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH]=eng_corr

              if Gname not in alldata[ads].keys():
                  alldata[ads][Gname]={}
              if toads not in alldata[ads][Gname].keys():
                  alldata[ads][Gname][toads]={}

              alldata[ads][Gname][toads][pH] = alldata[ads][Ename][toads][pH].copy()
              if 'ddag_IS' in Gname or 'ddag_FS' in Gname:
                  alldata[ads][Gname][toads][pH][1] += alldata[ads]['free_en_corr_%s'%facet]
              else:
                  alldata[ads][Gname][toads][pH][1] += alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH]

              if not ien:
               if 'E_ddag_vs_pot_quad_%s'%facet in alldata[ads]:
                  if 'G_ddag_vs_pot_quad_%s'%facet not in alldata[ads]:
                      alldata[ads]['G_ddag_vs_pot_quad_%s'%facet]={}
                  if toads in alldata[ads]['E_ddag_vs_pot_quad_%s'%facet]:
                      if toads not in alldata[ads]['G_ddag_vs_pot_quad_%s'%facet]:
                          alldata[ads]['G_ddag_vs_pot_quad_%s'%facet][toads]={}
                      if pH in alldata[ads]['E_ddag_vs_pot_quad_%s'%facet][toads]:
                          #print(ads,toads)
                          alldata[ads]['G_ddag_vs_pot_quad_%s'%facet][toads] = alldata[ads]['E_ddag_vs_pot_quad_%s'%facet][toads].copy()
                          alldata[ads]['G_ddag_vs_pot_quad_%s'%facet][toads][pH][2] += alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH]

    print('Vibrations could not be found for: ', no_freen)

def get_data(potentials,facets,include_barriers=True,basepath=os.getcwd(),cavity_file='cavity_pot_2.out',
             barrier_path='/Users/geokast/SelectCO2/barriers', outfile = 'results/parsed_data.pckl',plot=False,alldata={},
             thermo_path='dft_data/thermodynamics/',vibfile='dft_data/vibrations.pckl',
             backup_vibfile='libs/Cu_surface_sampling_constq.txt'):

    #If a single facet is given as a string
    if isinstance(facets,str):
        facets = [facets]
    for facet in facets:
        print('-'*15)
        print(f'Parsing {facet}')
        facpath=thermo_path
        not_enough_potentials=[]
        for ads in os.listdir(facpath):
            final_paths=[]
            adspath=os.path.join(facpath,ads)
            if os.path.isdir(adspath):
                for inads in os.listdir(adspath):
                    if inads.split('_')[0] == ads:
                        final_paths.append(os.path.join(adspath,inads))

            if not len(final_paths):  continue

            if ads not in alldata.keys():
                alldata[ads]={}
            adsdata={}
            for final_path in  final_paths:
                sitename=final_path.split('/')[-1]
                adsdata[sitename] = {}
                for adsfile in os.listdir(final_path):
                    for potential in np.around(potentials,2):
                        if len(adsfile.split('_')) > 1:
                            if '_'.join(adsfile.split('_')[:2]) == 'Pot_%1.2f'%potential:
                                if 'E_%s'%facet not in adsdata[sitename].keys():
                                    adsdata[sitename]['E_%s'%facet]={}
                                    adsdata[sitename]['E_C_%s'%facet]={}
                                if 'ne_%s'%facet not in adsdata[sitename].keys():
                                    adsdata[sitename]['ne_%s'%facet]={}
                                #print(sitename,facet,potential)
                                adsdata[sitename]['E_%s'%facet][potential] = read(final_path+'/'+adsfile).get_potential_energy()
                                adsdata[sitename]['ne_%s'%facet][potential] = read(final_path+'/'+adsfile).calc.results['ne']
                                elpot=read(final_path+'/'+adsfile).calc.results['electrode_potential']
                                adsdata[sitename]['E_C_%s'%facet][adsdata[sitename]['ne_%s'%facet][potential]] =\
                                        read(final_path+'/'+adsfile).get_potential_energy() -\
                                        adsdata[sitename]['ne_%s'%facet][potential]*elpot

                                if ads == 'clean':
                                    adsdata[sitename]['cell'] = read(final_path+'/'+adsfile).cell
                    if adsfile == cavity_file:
                        adsdata[sitename]['cavity_%s'%facet] = np.loadtxt(final_path+'/'+adsfile)

            # Find the most stable site from selection
            #if facet == '211' and ads == 'CO':
                #testprint(alldata[ads])
            #    alldata[ads].update(find_most_stable_site(adsdata,facet,potentials).copy())
            #    testprint(ads,alldata[ads].keys())
            #    das
            #else:
            alldata[ads].update(find_most_stable_site(adsdata,facet,
                potentials,plot=plot,ads=ads,not_enough_potentials=not_enough_potentials).copy())

        E_C_clean_fit_data=[]
        for ne in alldata['clean']['E_C_%s'%facet].keys():
            E_C_clean_fit_data.append([ne,alldata['clean']['E_C_%s'%facet][ne]])
            #print(ne,alldata['clean']['E_C_%s'%facet][ne])
        E_C_clean_fit_data=np.array(E_C_clean_fit_data)
        coeff,d=curve_fit(quad_fun,E_C_clean_fit_data[:,0],E_C_clean_fit_data[:,1])
        alldata['clean']['E_C_vs_ne_%s'%facet]=coeff
        alldata,not_enough_potentials=fit_potential_response(alldata,facet,potentials,include_barriers=False,
                plot=False,specific_ads='clean',quad_fit=True,not_enough_potentials=not_enough_potentials)

        #Thermodynamics
        for ads in alldata.keys():
            #Loop over lots of potentials for checking inhomogeneous potential grids
            for potential in np.around(potentials,2):
                if 'E_%s'%facet in alldata[ads].keys() and ads != 'clean':
                    #Reference energies to the clean slab
                    try:
                        #alldata[ads]['E_%s'%facet][potential] -= \#alldata['clean']['E_%s'%facet][potential]
                        alldata[ads]['E_%s'%facet][potential] -= \
                                alldata['clean']['E_vs_pot_%s'%facet][0]*potential**2+\
                                alldata['clean']['E_vs_pot_%s'%facet][1]*potential+\
                                alldata['clean']['E_vs_pot_%s'%facet][2]
                    except KeyError:
                        #print('Potential %1.2f seems  to be missing for adsorbate %s'%(potential,ads))
                        pass
                    else:
                        ads_short=ads.lstrip('md-').lstrip('bdo-').replace('2','H')
                        ads_short=ads_short.replace('3','HH').replace('-uw','')
                        ads_short=ads_short.replace('-hollow','').rstrip('old').lstrip('x-')
                        alldata[ads]['E_%s'%facet][potential] -= get_reference_energies(ads_short,code='GPAW')
                    #Add canonical energies
                    try:
                        alldata[ads]['E_C_%s'%facet][alldata[ads]['ne_%s'%facet][potential]] -=\
                                alldata['clean']['E_C_vs_ne_%s'%facet][0]*alldata[ads]['ne_%s'%facet][potential]**2+\
                                alldata['clean']['E_C_vs_ne_%s'%facet][1]*alldata[ads]['ne_%s'%facet][potential]+\
                                alldata['clean']['E_C_vs_ne_%s'%facet][2]
                    except KeyError:
                        pass
                    else:
                        alldata[ads]['E_C_%s'%facet][alldata[ads]['ne_%s'%facet][potential]] -=\
                                get_reference_energies(ads_short,code='GPAW')

        if include_barriers:
            alldata = get_barriers(alldata,barrier_path,facet,plot_charge_transfer=plot)
        add_gas_phase_data_to_dict(alldata,facet)

        alldata,not_enough_potentials=fit_potential_response(alldata,facet,potentials,include_barriers=True,plot=plot,not_enough_potentials=not_enough_potentials)
        if len(not_enough_potentials):
            print('More potential should be calculated for: ',not_enough_potentials)
        read_vibrational_frequencies(alldata,None,backup_vibfile,facet,vibfile=vibfile)
        add_vibrational_free_energy_corrections(alldata,facet)
        if plot:
            plot_Eddag_from_specific_IS(alldata,'G_ddag',facet,potentials,True,False)

        #testcoeff=alldata['CO']['E_ddag_vs_pot_%s'%facet]['COH']['acid']
        #print('getdata after addvibs CO COH',testcoeff,testcoeff[0]*4.4+testcoeff[1])

    import pickle
    out=open(outfile,'wb')
    pickle.dump(alldata,out)

    return alldata

def add_gas_phase_data_to_dict(alldata,facet=None,temperature=298.15,references={'C':'CO'}):
    #Reference gas phase molecules
    from ase.thermochemistry import  IdealGasThermo
    from ase.units import invcm
    gases=['H2_g','H2O_g','CO_g','CH4_g','C2H4_g','CH3CH2OH_g','CH3COOH_g','CH3CHO_g']
    for gas in gases:
        alldata[gas]={}

    alldata['H2_g']['E']=0
    alldata['H2O_g']['E']=0
    alldata['CO_g']['E']=0

    alldata['H2_g']['vibs']=[0, 123.3, 182.2, 304.6, 427.9, 4470.7]
    alldata['CO_g']['vibs']=[0, 0, 78.6, 263.9, 284.3, 2115.0]
    alldata['H2O_g']['vibs']=[0, 0, 121.9, 259.2, 290.5, 358.9, 1623.4, 3753.9,  3874.9]

    alldata['CH4_g']['E']=-2.486809
    alldata['CH4_g']['vibs']=[0, 0, 0,  0,  83.0, 143.7, 1288.3, 1302.9, 1304.3, 1504.7, 1515.7, 3022.0, 3108.2, 3121.0, 3126.1]
    alldata['C2H4_g']['E']=-2.76975805
    alldata['C2H4_g']['vibs']=[0, 0, 0, 0, 0, 0, 826.36096307,  959.45228605,  960.34934469, 1053.97242644, 1229.69837936, 1348.85228724, 1472.24547713, 1668.96347041, 3105.59055667, 3125.65012816, 3173.82607546, 3196.77938303]

    alldata['CH3CH2OH_g']['E']=-3.255991999999992
    alldata['CH3CH2OH_g']['vibs'] = [0,  0,   0,  134.7,  149.6,  181.7,  334.3,  384.9,  495.9,  832.0, 907.2,  992.8, 1081.6, 1157.1, 1271.7, 1295.2, 1403.1, 1426.8, 1465.6, 1502.8, 1519.1, 2959.8, 2971.0, 3032.1, 3073.8, 3090.5, 3799.5]

    alldata['CH3COOH_g']['E']=-2.7303311098416003
    alldata['CH3COOH_g']['vibs']=[0,0,0,0,92.5,254.2, 318.8,  413.9,525.4,623.6, 654.9, 867.1, 943.1,1050.2,1122.5, 1294.8, 1365.3, 1414.4, 1431.9, 1774.1, 3014.7, 3085.4, 3154.0, 3705.6]
    alldata['CH3CHO_g']['E']=-2.4735291705117106
    alldata['CH3CHO_g']['vibs']=[0,0,0,0,92.9,198.7,225.1,550.9,691.1,916.3,1074.4,1085.2,1332.4,1381.1,1424.6,1457.6,1768.3,2812.4,2985.8,3079.3,3103.4]

    alldata['CO_g'].update({'pressure':101325,'geometry':'linear','symmetry': 1})
    alldata['H2_g'].update({'pressure':101325,'geometry':'linear','symmetry': 1})
    alldata['H2O_g'].update({'pressure':0.035*101325,'geometry':'nonlinear','symmetry':2})

    alldata['C2H4_g'].update({'pressure': 1,'geometry':'nonlinear','symmetry': 4})
    alldata['CH3CH2OH_g'].update({'pressure': 1, 'geometry':'nonlinear','symmetry':1})
    alldata['CH4_g'].update({'pressure': 1,'geometry':'nonlinear','symmetry':12})
    alldata['CH3COOH_g'].update({'pressure': 1, 'geometry':'nonlinear','symmetry':1})
    alldata['CH3CHO_g'].update({'pressure': 1, 'geometry':'nonlinear','symmetry':1})

    for ads in gases:
        gibbs = IdealGasThermo(vib_energies = np.array(alldata[ads]['vibs'])*invcm,
                                    geometry=alldata[ads]['geometry'],
                                    spin=0,
                                    symmetrynumber=alldata[ads]['symmetry'],
                                    atoms=read('/Users/geokast/SelectCO2/endstates/gas_geometries/'+ads.rstrip('_g')+'.traj'))

        alldata[ads]['free_en_corr'] = gibbs.get_gibbs_energy(
                                                pressure=alldata[ads]['pressure'],
                                                temperature=temperature,
                                                verbose=False)

        alldata[ads]['G'] = alldata[ads]['E']+alldata[ads]['free_en_corr']


    vibnames={'CH4_g':'CHHH','C2H4_g':'CCHHHH','CH3CH2OH_g':'CHHHCHHOH','CH3COOH_g':'CHHHCOOH','CH3CHO_g':'CHHHCHO','H2_g':'HH'}
    for ads in ['CH4_g','C2H4_g','CH3CH2OH_g','CH3COOH_g','CH3CHO_g','H2_g']:
        alldata[ads]['G'] -= get_reference_vibrational_contribution(vibnames[ads],references=references)

    #print('Equilibrium potential of CH4:',-alldata['CH4_g']['G']/6)#*0.255)
    #print('Equilibrium potential of C2H4_g:',-alldata['C2H4_g']['G']/8)#+8*0.14)
    #print('Equilibrium potential of CH3CH2OH_g:',-alldata['CH3CH2OH_g']['G']/8)#+8*0.14)
#    das



def find_most_stable_site(data,facet,potentials,potential_to_check=2.9,plot=False,
        plotdir='results/E_on_varying_sites/',ads=None,quad_fit=True,not_enough_potentials=None):
    E_at_pots=[]
    if not_enough_potentials is None:
        not_enough_potentials=[]
    data,not_enough_potentials=fit_potential_response(data,facet,potentials,include_barriers=False,plot=False,quad_fit=quad_fit,not_enough_potentials=not_enough_potentials)
    coeffs=[]
    for site in data:
        if 'E_vs_pot_%s'%facet not in data[site].keys():continue
        coeff=data[site]['E_vs_pot_%s'%facet].copy()
        #the "fit_potential_response" function defaults to "E_vs_pot" for the fit
        #It's manually renamed here
        data[site]['E_abs_vs_pot_%s'%facet]=data[site]['E_vs_pot_%s'%facet].copy()
        del data[site]['E_vs_pot_%s'%facet]

        if quad_fit:
            E_at_pots.append([site,potential_to_check**2*coeff[0]+potential_to_check*coeff[1]+coeff[2]])
        else:
            E_at_pots.append([site,potential_to_check*coeff[0]+coeff[1]])
        coeffs.append([site,coeff])

    if not len(E_at_pots): return data[site]
    i_most_stable_site=np.argsort(np.array(E_at_pots)[:,1])[-1]

    if plot:
        colors=cm.gist_ncar(np.linspace(0,0.5,len(coeffs)))
        for ic,coeff in enumerate(coeffs):
            #print(coeff,data[coeff[0]]['E_%s'%facet])
            E_at_site=[]
            for pot in np.linspace(2.4,3.4,5):
                E_at_site.append([pot,coeff[1][0]*pot**2+coeff[1][1]*pot+coeff[1][2]])
            E_at_site=np.array(E_at_site)
            if len(coeff[0].split('_')) > 3:
                label = ' '.join(coeff[0].split('_')[3:])
            else:
                label = ' '.join(coeff[0].split('_')[1:])
            plt.plot(E_at_site[:,0],E_at_site[:,1],'-',
                         label=label,
                         color=colors[ic])
            points=[]
            for pot in data[coeff[0]]['E_%s'%facet]:
                points.append([float(pot),data[coeff[0]]['E_%s'%facet][pot]])
            points=np.array(points)
            plt.plot(points[:,0],points[:,1],'o',color=colors[ic],markeredgecolor='k')
        plt.ylabel('$\Omega$ [eV]')
        plt.xlabel('Work function [eV]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(plotdir+'E_vs_site_%s_%s.pdf'%(ads,facet))
        plt.close()

    return data[E_at_pots[i_most_stable_site][0]]


def get_barriers(alldata,barrier_path,facet,plot_charge_transfer=False):
    barrier_dirs=[]
    #Check which barriers might be there
    for bardir in os.listdir(barrier_path):
        if (len(bardir.split('_')) > 2 and bardir != 'CO_to_OCCO'):
            if bardir.split('_')[1] == 'to':
                if facet in os.listdir(barrier_path+'/'+bardir):
                    barrier_dirs.append(barrier_path+'/'+bardir+'/'+facet)

    #Check at which potential the barriers have been calculated
    for bardir_iter in barrier_dirs:
      #Identify IS and FS
      ads=bardir_iter.split('/')[-2].split('_')[0]
      toads=bardir_iter.split('/')[-2].split('_')[2]

      if ads not in alldata.keys():
            if 'bdo-'+ads in alldata.keys():
                ads='bdo-'+ads
            elif 'md-'+ads in alldata.keys():
                ads='md-'+ads
            elif ads == 'COCOH':
                ads='COH'
            elif ads == 'COCHO':
                ads='CHO'
            else:
                print('Could not find the initial state for barrier %s'%bardir_iter.split('/')[-2])
                continue


      toads_short=toads.lstrip('md-').lstrip('bdo-')
      ads_short=ads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')
      toads_short2=toads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')
      letter_diff=toads_short2
      for letter in ads_short:
          letter_diff=letter_diff.replace(letter,'',1)

      for outname in ['E_ddag_%s'%facet,'E_ddag_IS_%s'%facet,'E_ddag_FS_%s'%facet]:
          if outname not in alldata[ads].keys():
              alldata[ads][outname]={}
          alldata[ads][outname][toads_short]={}

      # Determine all directories barriers have been calculated in for the pair
      # of IS and FS
      bardirs=[]
      for pH in ['base','acid','chemical']:
       charges,energies={},{}
       if pH == 'acid':
           if 'acidic' not in os.listdir(bardir_iter):  continue
           bardirs.append(bardir_iter+'/acidic')
       elif pH == 'chemical':
           if 'chemical' not in os.listdir(bardir_iter):  continue
           bardirs.append(bardir_iter+'/chemical')
       else:
           bardirs.append(bardir_iter)

       if ads_short == 'OHCCHH' and toads_short == 'OCHCH3':
           bardirs[-1]+='/desorption/'
       #elif all([ads_short == 'H',toads_short == 'H2',facet=='100',pH=='base']):
       #    bardirs.append(bardir_iter)
       #    bardirs[-2]+='/H-down/'
          # bardirs[-1]+='/H-up/'


      for bardir in bardirs:
       print(ads_short,toads_short,bardir)

       #Read pH from path
       if 'acidic' in bardir:
           pH='acid'
       elif 'chemical' in bardir:
           pH='chemical'
       else:
           pH='base'

       # Check which potentials have been calculated
       potentials=[]
       for potdir in os.listdir(bardir):
            if len(potdir.split('_')) > 1:
                if potdir.split('_')[0] == 'pot':
                    if 'neb_GC_final_climbed.traj' in os.listdir(bardir+'/'+potdir):
                        potentials.append(potdir.split('_')[1])


       #Determine the difference of IS and FS for formation energies
       for potential in sorted(potentials):
          #print(bardir)
          atoms_list=read(bardir+'/pot_%s/neb_GC_final_climbed.traj@:'%potential)

          for state in ['barrier','IS','FS']:
              outname='E_pFS_%s'%facet
              if state == 'barrier':
                  Eddag=max([atoms.get_potential_energy() for atoms in atoms_list])
                  outname='E_ddag_%s'%facet
              elif state == 'IS':
                  Eddag=atoms_list[0].get_potential_energy()
                  outname='E_ddag_IS_%s'%facet
              elif state == 'FS':
                  Eddag=atoms_list[-1].get_potential_energy()
                  outname='E_ddag_FS_%s'%facet
    #          if pH == 'acid':
    #              print('1',ads,potential,Eddag,alldata['clean']['E_%s'%facet][float(potential)])
              if float(potential) in alldata['clean']['E_%s'%facet]:
                  Eddag-=alldata['clean']['E_%s'%facet][float(potential)]
              elif 'E_vs_pot_%s'%facet in alldata['clean']:
                  Eddag-= alldata['clean']['E_vs_pot_%s'%facet][0]*float(potential)**2+\
                          alldata['clean']['E_vs_pot_%s'%facet][1]*float(potential)+\
                          alldata['clean']['E_vs_pot_%s'%facet][2]
              else:
                  print('Potential %s of the clean slab could not be found for barrier from %s to %s.'%(potential,ads,toads))
                  continue

              #It's always the IS we reference to in alkaline
              #TODO: Check whether this works

              #Protonation
              if letter_diff == 'H':
              #    if any(['acidic' in bardir, 'chemical' in bardir]):
              #     reference_string=toads_short2#+'HHO'
              #    else:
                  if pH == 'base':
                   reference_string=ads_short#+'HHO'
                  elif pH in ['acid','chemical']:
                   reference_string=toads_short2#+'HHO'

              else:
    #              if 'chemical' in bardir:
    #                  raise NotImplementedError('Chemical OH desorption is not '
    #                          'implemented')
    #              elif 'acidic' in bardir:
    #                reference_string=toads_short+'HHO'
    #              else:
    #                reference_string=ads_short
                  if pH == 'base':
                    reference_string=ads_short
                  elif pH in ['acid']:
                    reference_string=toads_short+'HHO'
                  elif pH in ['chemical']:
                      raise NotImplementedError('Chemical OH desorption is not '
                              'implemented')

              if toads_short in ['OCCO']:
                  Eddag -= get_reference_energies(toads_short,code='GPAW')
              elif (toads_short in ['OCCOH'] and ads == 'COH') or\
                      (toads_short in ['OCCHO'] and ads == 'CHO'):
                  Eddag -= get_reference_energies(toads_short,code='GPAW')
              elif reference_string == 'clean':
                  pass
              else: # len(ads_short) > len(toads_short):
                  Eddag -= get_reference_energies(reference_string,code='GPAW')

              if state == 'barrier':
                  try:
                      charges[potential]= [i.calc.results['ne'] for i in atoms_list]
                      energies[potential] = [i.calc.results['energy'] for i in atoms_list]
                  except:
                      pass

              if pH not in alldata[ads][outname][toads_short].keys():
                  alldata[ads][outname][toads_short][pH] = {}
              alldata[ads][outname][toads_short][pH][potential]=Eddag
          #if ads=='COCO':
          #  if 'Eddag' not in alldata['CO']:
          #      alldata['CO']['E_ddag_%s'%facet]={toads: {pH: {potential: Eddag}}}
          #  alldata['CO']['E_ddag_%s'%facet][toads][pH][potential]=Eddag
#       if facet == '211':
#           print(facet,ads,toads,potential,alldata[ads]['E_ddag_%s'%facet][toads])
          if len(charges.keys()) and plot_charge_transfer:
              plot_ne_over_band(charges,energies,ads,toads_short,bardir_iter,pH)
      #print(bardir,potentials)
    #das
    #print(alldata['clean'])
    return alldata

def plot_ne_over_band(charges,energies,ads,toads,bardir,pH):
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('#image')
    ax1.set_ylabel(r'$\Delta$E [eV]')
    ax2=ax1.twinx()
    ax2.set_ylabel('Charge transfer [e]')
    for pot in charges.keys():
        ax2.plot(np.array(charges[pot])-charges[pot][0],'--')
        ax1.plot(np.array(energies[pot])-energies[pot][0],label='WF=%1.2feV'%float(pot))
    ax1.set_title(bardir.split('/')[-1])
    fig.tight_layout()
    fig.legend()
    fig.savefig('results/charge_transfer_along_band/charge_transfer_'+bardir.split('/')[-1]+'_'+pH+'.pdf')
    plt.close()

def read_vibrational_frequencies(alldata,line,backup_vibfile=None,facet='100',vibfile='dft_data/vibrations.pckl',use_HHH_vibs=False,substrates=['Cu'],no_water_layer=False):
    # "no_water_layer" has been added for the calculation without water layer in the barrier calculations

    try:
        allvibs=pickle.load(open(vibfile,'rb'))
    except FileNotFoundError:
        print('WARNING! Could not load the vibrations pickle file. Free energies not accessible')
        return

    #print(allvibs.keys())
    #print(alldata.keys())

    missing_vibs=[]
    used_vibs=[]
    for ads_long in alldata.keys():
     if no_water_layer and ads_long == 'clean_slab': continue
     for pH in ['base']:
        no_vibs=True
        if ads_long[-2:] == '_g': continue
        ads=ads_long.lstrip('md-').lstrip('bdo-')
        # Read vibrational frequencies of thermodynamics
        if ads_long in allvibs.keys() and not use_HHH_vibs:
         if facet in allvibs[ads_long].keys():
          if 'vibs' in allvibs[ads_long][facet].keys():
           if pH in allvibs[ads_long][facet]['vibs'].keys():
            if 'vibs_%s'%facet not in alldata[ads_long].keys():
                alldata[ads_long]['vibs_%s'%facet] = {}
            alldata[ads_long]['vibs_%s'%facet] = [i/units.invcm for i in allvibs[ads_long][facet]['vibs'][pH]]
            no_vibs=False
        #    continue

        # If the vibs are found everything is good and we go to the next
#        print(ads_long)
        if not no_vibs:
            used_vibs.append(ads_long)
            continue

        if backup_vibfile:
            missing_vibs.append(ads_long)
            #print(ads_long+' vibs have not been found in my vibs')
            viblines=open(backup_vibfile,'r').readlines()[1:]
            for ivibline in viblines:
                if ivibline.split()[0] in substrates+['None']:
                    #if ivibline.split()[1] in [facet,'gas']:
                    if ivibline.split()[1] in [facet]:
                        if ivibline.split()[2] == ads:
                           vibline=ivibline
                           break
            else:
      #          print(ads + ' is not in any vibfile')
                continue
                #ads_and_electron[ads]['vibs_%s'%facet] = []
#                return []
        else:
            vibline=line

        if vibline is None:
            continue
        freq_inline=[None,None]
     #   print(ads_long)
        for isplit,splitline in  enumerate(vibline.split()):
            if splitline[0] == '[':
                freq_inline[0]=isplit
            elif splitline[-1] == ']':
                freq_inline[1]=isplit+1
                break

        if None not in freq_inline:
            frequencies = [float(vib.replace(',','').replace('[','').replace(']',''))
                    for vib in vibline.split()[freq_inline[0]:freq_inline[1]]]

        else:
            print('No frequencies given for '+ads)
            frequencies=[]

        alldata[ads_long]['vibs_%s'%facet] = frequencies

    print(f'Missing vibrations for {facet}:', missing_vibs)

    if no_water_layer:
        bars=[barads
              for barads in allvibs.keys()
              if '-' in barads]
    else:
        bars=[barads
              for barads in allvibs.keys()
              if '_to_'  in barads]

    for bar in bars:
        for pH in ['base','acid','chemical']:
            if facet not in allvibs[bar].keys(): continue
            if 'vibs' not in allvibs[bar][facet].keys(): continue
            if pH not in allvibs[bar][facet]['vibs'].keys(): continue
            if not len(allvibs[bar][facet]['vibs'][pH]): continue

            # The following "if no_water_layer" has been added for the calculation without water layer
            if no_water_layer:
              IS_FS={'OC-CO': ['COCO','OCCO'],
                    'OCCO-H': ['OCCO','OCCOH'],
                    'H-HCCO': ['HCCO','H2CCO'],
                    'HCCO-H': ['HCCO','HCCOH'],
                    'HOCCO-H':['OCCOH','HOCCOH'],
                    'OCC-OH':['OCCOH','CCO'],
                    'CO$_{2(g)}$_to_HCOO$^-_{(aq)}':['CO$_{2(g)}$','HCOO$^-_{(aq)}$'],
                    'OC-O':['CO$_{2(g)}$','CO2']}
              barads=IS_FS[bar][0]
              bartoads=IS_FS[bar][1]
            else:
                barads=bar.split('_')[0]
                bartoads=bar.split('_')[-1]
                bartoads=bartoads.lstrip('md-').lstrip('bdo-')
                if barads == 'COCOH':
                    barads='COH'
                elif barads == 'COCHO':
                    barads='CHO'
            #print(alldata.keys())
            if 'vibs_ddag_%s'%facet not in alldata[barads].keys():
                alldata[barads]['vibs_ddag_%s'%facet]={}

            if bartoads not in alldata[barads]['vibs_ddag_%s'%facet].keys():
                alldata[barads]['vibs_ddag_%s'%facet][bartoads]={}
            alldata[barads]['vibs_ddag_%s'%facet][bartoads][pH]=\
                    [i/units.invcm for i in allvibs[bar][facet]['vibs'][pH]]
            used_vibs.append(bar)

    unused_vibs=[]
    for vib in allvibs:
        if vib not in used_vibs and facet in allvibs[vib]:
            unused_vibs.append(vib)
    print(f'Unused vibrations for {facet}:', unused_vibs)

def plot_imaginary_frequencies(alldata, facet,outfilebase='results/img_frequencies.pdf',
        pHs=['base','acid','chemical']):

    for pH in pHs:
      outfile=outfilebase.rstrip('.pdf')+'_%s.pdf'%pH
      xlabels=[]
      for iads,ads in enumerate(alldata.keys()):
        if 'vibs_ddag_%s'%facet not in alldata[ads].keys(): continue
        for itoads,toads in enumerate(alldata[ads]['vibs_ddag_%s'%facet].keys()):
          #for pH in alldata[ads]['vibs_ddag_%s'%facet][toads].keys():
            if pH not in alldata[ads]['vibs_ddag_%s'%facet][toads]: continue
            if not np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads][pH][0]): continue

            plt.plot(len(xlabels),np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads][pH][0])
                    ,'o',color=get_intcolors(ads,toads),markeredgecolor='k',markersize=8)

            ads_short = ads.lstrip('md-').lstrip('bdo-')
            toads_short = toads.lstrip('md-').lstrip('bdo-')
            xlabels.append(get_intcolors(ads,toads,return_name=True))

      plt.plot(np.nan,np.nan,'o',color=get_intcolors('CO','CHO'),label='C protonation')
      plt.plot(np.nan,np.nan,'o'+get_intcolors('HOCCOH','CCOH'),label='OH desorption')
      plt.plot(np.nan,np.nan,'o'+get_intcolors('OCCO','OCCOH'),label='O protonation')
      plt.plot(np.nan,np.nan,'o'+get_intcolors('HCCO','H2CCO'),label='Anionic desorption')
      plt.plot(np.nan,np.nan,'o'+get_intcolors('clean','H'),label='Other')
      #plt.legend(bbox_to_anchor=(1,1))
      plt.legend()
      plt.xticks(np.arange(len(xlabels)),xlabels,rotation='vertical')
      plt.ylim=[0,1800]
      plt.ylabel('Imaginary frequency [cm$^{-1}$]',fontsize=14)
      plt.tight_layout()
      #plt.show()
      plt.savefig(outfile)
      plt.close()

def apply_Wigner_correction(alldata,facets):
    from ase.units import invcm
    kb=8.617333262145e-5
    T=300.
    hbar=6.582119569e-16
    fig,ax=plt.subplots(1,2,sharex=True,figsize=(12,6))
    fig2,ax2=plt.subplots(1,2,sharex=True,figsize=(12,6))
    second_order_wigfacs,first_order_wigfacs,imgvibs=[],[],[]
    for facet in facets:
     for iads,ads in enumerate(alldata.keys()):
        if 'vibs_ddag_%s'%facet not in alldata[ads].keys(): continue
        for itoads,toads in enumerate(alldata[ads]['vibs_ddag_%s'%facet].keys()):
            if not np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads]['base'][0]): continue
            imgvib=np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads]['base'][0])
            print(ads,toads,imgvib)
            beta_vib=imgvib*invcm/(kb*T)
            print(imgvib,imgvib*invcm,beta_vib,imgvib*invcm/(2*np.pi*kb))

            #Tunneling factor from
            #https://pubs.rsc.org/en/content/articlelanding/2014/CP/C4CP03235G#!divAbstract
            #qi0=1/
            #vibs=np.array(alldata[ads]['vibs_ddag_%s'%facet][toads])


            if 0:
               for N in range(10,100):
                    A_N=1
                    for j in range(3,N-1):

                        eta_0_j = np.sqrt((2*j*np.pi*kb*T/hbar)**2+(imgvib*invcm/(hbar*2*np.pi))**2)
                        print(eta_0_j,kb*T/(hbar*2*np.pi*eta_0_j),hbar*2*np.pi*eta_0_j)
                        A_N*=kb*T/(hbar*2*np.pi*eta_0_j)
               print(A_N)
            #Rate from New Journal ofPhysics 12 (2010) 055002
            om_0=3000
            om_b=1000
            E_b=.7
            T1=300

            for om_b in np.linspace(1,2000,30):
            #for T1 in [100,300,1000]:
                k=om_b*invcm/(4*np.pi*hbar)*\
                  np.sinh(om_0*invcm/(2*kb*T1))/np.sin(om_b*invcm/(2*kb*T1))*np.exp(-E_b/(kb*T1))
                kTST=kb*T1/(2*np.pi*hbar)*np.exp(-E_b/(kb*T1))
                print(T1,k,kTST,k/kTST)
                plt.plot(om_b,k/kTST,'o')
            #plt.ylim([-10,10])
            plt.show()
            asd

            #for k in range(10):
            #    eta_0_k = np.sqrt((2*k*np.pi*kb*T/hbar)**2+(imgvib*invcm/(hbar*2*np.pi))**2)

            #wigfac = beta_vib/2/(np.sin(beta_vib/2))
            crossover_T=imgvib*invcm/(2*np.pi*kb)
            print('Crossover temperature: %s-%s%s, '%(ads,toads,facet), crossover_T)
            second_order_wigfac=1+1/24*(beta_vib)**2
            second_order_wigfacs.append(second_order_wigfac)
            first_order_wigfac=(beta_vib)/2* 1/(np.sin(beta_vib/2.))
            first_order_wigfacs.append(first_order_wigfac)
            imgvibs.append(imgvib)
            ax[1].plot(imgvib,second_order_wigfac,'o',color=get_intcolors(ads,toads))#,label='Second order Wigner')
            ax2[1].plot(imgvib,second_order_wigfac,'o',color=get_intcolors(ads,toads))#,label='Second order Wigner')
            ax[1].plot(imgvib,first_order_wigfac,'d',color=get_intcolors(ads,toads))#,label='Standard Wigner')
            ax2[0].plot(imgvib,first_order_wigfac,'d',color=get_intcolors(ads,toads))#,label='Standard Wigner')
            ax[0].plot(imgvib,crossover_T,'o',color=get_intcolors(ads,toads))#,label='Second order Wigner')
    imgvibs=np.array(imgvibs)
    sort_img=np.argsort(imgvibs)
    imgvibs=imgvibs[sort_img]
    first_order_wigfacs=np.array(first_order_wigfacs)[sort_img]
    second_order_wigfacs=np.array(second_order_wigfacs)[sort_img]

    #Imaginary vibration of crossover at given T
    crossover_vib=(2*np.pi*kb*T)/invcm
    ax[1].axvline(x=crossover_vib,linestyle='--',color='k')
    ax2[0].axvline(x=crossover_vib,linestyle='--',color='k')
    ax2[1].axvline(x=crossover_vib,linestyle='--',color='k')
    ax2[0].annotate(r'$\omega_c$',(1310,1),fontsize=20)
    ax2[1].annotate(r'$\omega_c$',(1310,1),fontsize=20)
#    ax[1].plot(imgvibs,first_order_wigfacs,'--k')
#    ax[1].plot(imgvibs,second_order_wigfacs,'-k')
    ax[0].axhline(y=T,linestyle='--',color='k')
    ax[1].annotate('T=%iK'%int(T),(1300,0.1),fontsize=20)
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','CHO'),label='C protonation')
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','COH'),label='O protonation')
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('COH','C'),label='OH desorption')
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('clean','H'),label='Other')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','CHO'),label='C protonation')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','COH'),label='O protonation')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('COH','C'),label='OH desorption')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('clean','H'),label='Other')
    ax[1].plot(np.nan,np.nan,'ok',label='Second order Wigner')
    ax[1].plot(np.nan,np.nan,'dk',label='Standard Wigner')
    ax2[1].set_title('Second order Wigner')
    ax2[0].set_title('Standard Wigner')
    ax[0].set_ylabel('T$_{cross} [K]$')
    ax[1].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax[0].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax2[1].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax2[0].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax[1].set_ylabel('$\kappa_t$')
    ax2[0].set_ylabel('$\kappa_t$')
    ax[1].set_ylim([0,5])

    ax[0].legend()
    ax2[0].legend()
    ax[1].legend()
    fig.tight_layout()
    fig2.tight_layout()
    fig.savefig('results/Wigner_CrossoverT_and_tunneling.pdf')
    fig2.savefig('results/Wigner_tunneling_first_and_second_order.pdf')
    #plt.show()

def plot_Eddag_from_specific_IS(alldata,enkey,facet,potentials,deeperkey=None,quad_fit=False,plot_separate=True,ylabel='Formation energy [eV]'):

    symbols=['+','s','d']
    enkey2=enkey+'_vs_pot'
    for iads,ads in enumerate(alldata.keys()):
     plotted=False
     for quad_fit in [False,True]:
         enkey2=enkey+'_vs_pot'
         if quad_fit: enkey2+='_quad'
         if ads[-2:] == '_g': continue
         if enkey2+'_%s'%facet not in alldata[ads].keys(): continue
#             ads in ['clean']): continue
         #if enkey+'_vs_pot_%s'%facet in alldata[ads].keys() and ads not in ['clean']:
         if deeperkey:
         #  print(ads,enkey)
           counter=0
           for itoads,toads in enumerate(alldata[ads][enkey2+'_%s'%facet]):
              if enkey2+'_%s'%facet in alldata[ads].keys():
                 counter+=1
           colors=cm.jet(np.linspace(0,1,counter+1))

           counter2=0
           for toads in alldata[ads][enkey2+'_%s'%facet].keys():
            for pH in alldata[ads][enkey2+'_%s'%facet][toads].keys():

             coeff=alldata[ads][enkey2+'_%s'%facet][toads][pH]

             fit = []
             #print(ads,toads,enkey2,coeff,quad_fit)
             if quad_fit:
                 for pot in np.linspace(potentials[0],potentials[-1],10):
                     fit.append([pot,np.array([pot**2,pot,1])@coeff])
             else:
                 for pot in np.linspace(potentials[0],potentials[-1],2):
                     fit.append([pot,np.array([pot,1])@coeff])
             fit=np.array(fit)

             if quad_fit:
                label=ads+'-'+toads+', '+pH+r',d%s/d$\Phi$=(%1.2f,%1.2f)'%(enkey[0],2*coeff[0],coeff[1])
                plt.plot(fit[:,0],fit[:,1],'--',color=colors[counter2%len(colors)],
                     label=label)
             else:
                label=ads+'-'+toads+', '+pH+r',d%s/d$\Phi$=%1.2f'%(enkey[0],coeff[0])
                plt.plot(fit[:,0],fit[:,1],'-',color=colors[counter2%len(colors)],
                     label=label)
             try:
                E_v_pot=get_E_vs_pot(alldata,ads,potentials,enkey+'_%s'%facet,[toads,pH])
             except:
                 pass
             else:
                if len(E_v_pot):
                 plt.plot(E_v_pot[:,0],E_v_pot[:,1],
                         symbols[counter2%len(symbols)],
                         #label=ads+'-'+toads+', '+pH+r',dE/d$\Phi$=%1.2f'%coeff[0],
                         color=colors[counter2%len(colors)])
             plotted=True
             counter2+=1
     if not plotted: continue
     plt.xlabel('Work function [eV]')
     plt.ylabel(ylabel)
     plt.title('Barriers from %s on facet %s'%(ads,facet))
     #plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
     #           mode="expand", borderaxespad=0, ncol=4,fontsize=3)
     plt.legend()
     plt.tight_layout()
     plt.savefig('results/E_ddag_vs_pot_from_specific_ads/%s_vs_pot_%s_%s.pdf'%(enkey,ads,facet))
     plt.close()

def testprint(*string):
    print(*string)
