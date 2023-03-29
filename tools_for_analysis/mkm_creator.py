import numpy as np
import math
import sys,os
import pickle as pckl
from matplotlib import pyplot
from tools_for_analysis.sjm_analyse_tools import get_unique_TS_name,get_intcolors
#from sjm_analyse_tools import get_unique_TS_name,get_intcolors

def create_catmap_mkm(ads_and_electron,facet,dqdphi=20,pzc=-0.5,
        base_infile='/Users/geokast/SelectCO2/endstates/catmap/mkm_files/base_mkm.mkm',barriers=None,outfile=None):

    if outfile is None:
        outfile = '/Users/geokast/SelectCO2/endstates/catmap/mkm_files/%s.mkm'%facet

    inlines=open(base_infile,'r').readlines()

    for iline,line in enumerate(inlines):
        if line.split()[1:3] == ['Reaction','network']:
            rxn_line = iline+1
            break

    cat_spd_string=_create_species_definitions(ads_and_electron,facet,dqdphi,pzc)
    cat_rxn_string=_create_rxn_expressions(ads_and_electron,facet,barriers)

    out=open(outfile,'w')
    out.write(''.join(inlines[:rxn_line]))
    out.write(cat_rxn_string)
    out.write(''.join(inlines[rxn_line:]))
    out.write(cat_spd_string)
    out.write("\ndata_file = '%s.pkl'\n"%facet)
    out.close()
    print('-'*13)
    print('\nCatmap mkm file written in %s\n'%outfile)

    #_write_mkm_file()

def _add_sigma_params(ads_and_electron,facet,cat_string,dqdphi,pzc):

    cat_string+="\nsigma_input = ['CH', %i]\nUpzc = %1.2f\n"%(dqdphi,pzc)
    for ads in ads_and_electron.keys():
        #print(ads_and_electron[ads])
        if 'E_%s'%facet in ads_and_electron[ads].keys() and\
            'E_vs_pot_%s'%facet in ads_and_electron[ads].keys() and\
                ads_and_electron[ads]['n_sites'] > 0:
            #cat_string+="species_definitions['%s_%s']['sigma_params']=%s\n"%(ads,facet,ads_and_electron[ads]['dedphi_%s'%facet],ads_and_electron[ads]['E_%s'%facet])
            cat_string+="species_definitions['%s_%s']['sigma_params']=[%s]\n"%(ads,facet,','.join([str(i) for i in ads_and_electron[ads]['E_vs_pot_%s'%facet]]))
    return cat_string



def _create_species_definitions(ads_and_electron,facet,dqdphi,pzc):
    written_barriers=[]
    #Collect mono and bidentate adsorbates
    md,bd=[],[]
    cat_string=""
    cat_string=add_base_species_definitions(cat_string,facet)
    for ads_in in ads_and_electron.keys():
        ads=ads_in.replace('bdo-','').replace('md-','')
        if ads_and_electron[ads_in]['n_sites']==1:
            md.append(ads+'_'+facet)
        elif ads_and_electron[ads_in]['n_sites']==2:
            bd.append(ads+'_'+facet)
        else:
            print('Weird! %s has a unsuspected dentation'%ads_in)

        if 'E_ddag_vs_pot_%s'%facet in ads_and_electron[ads_in]:
            for toads in ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet].keys():
              if toads  in ads_and_electron.keys():
                #print(TS)
                #toads=TS.split('-')[1]
                TSname,written_barriers,dummy2 = get_unique_TS_name(ads,toads,'base',written_barriers)
                if any(np.array([ads_and_electron[ads_in]['n_sites'],ads_and_electron[toads]['n_sites']]) == 2):
                    bd.append(TSname+'_'+facet)
                else:
                    md.append(TSname+'_'+facet)


    #Monodentates treated like CO
    cat_string+="\n\nfor sp in ["
    for ads in md: cat_string+="'%s',"%ads
    cat_string+="]:\n"
    cat_string+="    species_definitions[sp] = species_definitions['CO_%s'].copy()\n"%facet

    #Bidentates treated like OCCO
    cat_string+="\nfor sp in ["
    for ads in bd: cat_string+="'%s',"%ads
    cat_string+="]:\n"
    cat_string+="    species_definitions[sp] = species_definitions['OCCO_%s'].copy()\n"%facet

    cat_string = _add_sigma_params(ads_and_electron,facet,cat_string,dqdphi,pzc)

    return cat_string


def _create_rxn_expressions(ads_and_electron,facet,barriers):
    enstring='E_%s'%facet
    print('-'*13)
    print('Creating catmap reaction expressions')
    #'H2O_g + ele_g + *_t  <-> H*_t + OH_g; beta=0.7',	# Volmer 1
    catmap_string=''
    counter=0
    written_barriers=[]
    for ads_in in ads_and_electron.keys():
        for toads_in in ads_and_electron[ads_in]['leads_to']:
            if enstring not in ads_and_electron[ads_in].keys():
                print('Reactant %s doesnt seem  to have a formation energy'%ads_in)
                continue
            if toads_in not in ads_and_electron.keys():
                print('Product %s is not in the given adsorbate list'%toads_in)
                continue
            if enstring not in ads_and_electron[toads_in].keys():
                print('Product %s doesnt seem to have a formation energy'%toads_in)
                continue
            ads=ads_in.replace('bdo-','',1).replace('md-','',1)
            toads=toads_in.replace('bdo-','',1).replace('md-','',1)


            elem_ads={}
            elem_toads={}

            #Make sure gas_phase species are not bound
            site={}
            for spec in [ads_in,toads_in]:
                if ads_and_electron[spec]['n_sites']:
                    site[spec]=facet
                else:
                    site[spec]='g'


            # Define the catmap mkm lines
            catmap_line,written_barriers,counter=\
                    get_catmap_base_line(ads_and_electron,ads_in,ads,toads,
                            toads_in,facet,barriers,site,counter,
                            written_barriers)
            catmap_string+=catmap_line

    catmap_string=finalize_rxn_expression(catmap_string,facet)
    return catmap_string

def get_catmap_base_line(ads_and_electron,ads_in,ads,toads,toads_in,facet,barriers,site,counter,written_barriers):
    catmap_line=''#string+=ads+'_'+facet
    if barriers:
        OH_des = {'base':'ele_g <-> %s <-> %s_%s + OH_g; beta=%s ',
                'acid': 'H_g + ele_g <-> %s <-> %s_%s + H2O_g;beta=%s'}
        H_ads = {'base':'ele_g + H2O_g <-> %s <-> %s_%s + OH_g; beta=%s',
                'acid': 'ele_g + H_g <-> %s <-> %s_%s + H2O_g; beta=%s'}
         #barr = ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet][toads]
    else:
        #OH_des = {'base': 'ele_g <-> %s_%s + OH_g'%(facet,toads,site[toads_in]),
        #        'acid': 'H_g + ele_g <-> %s_%s + H2O_g'%(facet,toads,site[toads_in])}
        #H_ads = {'base':'ele_g + H2O_g <-> %s_%s + OH_g'%(facet,toads,site[toads_in]),
        #         'acid':'ele_g + H_g <-> %s_%s + H2O_g'%(facet,toads,site[toads_in])}
        OH_des = {'base': 'ele_g <-> %s_%s + OH_g'%(toads,site[toads_in]),
                'acid': 'H_g + ele_g <-> %s_%s + H2O_g'%(toads,site[toads_in])}
        H_ads = {'base':'ele_g + H2O_g <-> %s_%s + OH_g'%(toads,site[toads_in]),
                 'acid':'ele_g + H_g <-> %s_%s + H2O_g'%(toads,site[toads_in])}

    if barriers:
     #Check if barriers have beeen calculated
     TStype=get_intcolors(ads,toads,return_type=True)
     for pH in ['base','acid']:
      if 'E_ddag_vs_pot_%s'%facet in ads_and_electron[ads_in].keys():
           if toads_in in ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet].keys():
             # XXX: Assuming linear fit of beta
             if pH in ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet][toads_in]:
                 beta=ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet][toads_in][pH][-2]
             else:
                 continue
           else:
           #    print(ads,ads_and_electron[ads_in])
               print('Beta for %s to %s not found'%(ads_in,toads_in))
               beta=0.5
           TSname,written_barriers,namediff=get_unique_TS_name(ads_in,toads_in,pH,written_barriers)
           TSname+='_%s'%facet
           if TStype in ['Cprot','Oprot','else']: catmap_line += H_ads[pH]%(TSname,toads,site[toads_in],beta)
#           elif TStype in ['Oprot']: catmap_line += O_H_ads%(TSname,toads,site[toads_in],beta)
           elif TStype in ['OHdes']: catmap_line += OH_des[pH]%(TSname,toads,site[toads_in],beta)

           #OH_des = OH_des%(get_unique_TS_name(ads_in,toads_in,pH,written_barriers),toads,site[toads_in],beta)
           #O_H_ads = O_H_ads%(ads,facet,toads,site[toads_in],beta)

     #If not calculated use assumed barriers
      elif not isinstance(barriers,dict):
         print('Barriers have to be given as a dict! Aborting')
         sys.exit()
      else:
         TSname='^%seV_%s'
         if TStype in ['Cprot','Oprot','else']:
             catmap_line += H_ads[pH]%(TSname%(barriers['C-H'],facet),toads,site[toads_in],0.5)
#         elif TStype in ['Oprot']:
#            catmap_line += 'ele_g + H2O_g <-> ^%1.4feV_%s <-> %s_%s + OH_g'%(barriers['CO-H'],facet,toads,site[toads_in])
         elif TStype in ['OHdes']:
             catmap_line += OH_des[pH]%(TSname%(barriers['COH-H'],facet),toads,site[toads_in],0.5)
#         elif TStype in ['else']:
#             catmap_line += 'ele_g + H2O_g <-> ^%1.4feV_%s <-> %s_%s + OH_g'%(barriers['C-H'],facet,toads,site[toads_in])


      #print(ads,toads,TStype)
      #   catmap_line+=C_H_ads
      #elif TStype == 'Oprot':
      #   catmap_line+=O_H_ads
      #elif TStype == 'OH_des':
      #   catmap_line+=OH_des

      if len(catmap_line) > 0:
                if ads == 'clean':
                    catmap_line2='*_%s + '%(site[toads_in])+catmap_line
                else:
                    catmap_line2='%s_%s + '%(ads,site[ads_in])+catmap_line
                counter+=1
                return "\'%s\',   #%i\n"%(catmap_line2,counter),written_barriers,counter
      else:
                print('Check:', ads,toads)
                return '',written_barriers,counter
     #return catmap_string


            #XXX: Make this more general


def finalize_rxn_expression(catmap_string,facet):
    catmap_string+="\
'H3CCH2O_{fac} + ele_g + H2O_g <-> ^0.5000eV_{fac} <-> CH3CH2OH_g + OH_g + 2*_{fac}',   #21\n\
'OH_100 + ele_g <-> ^0.2000eV_100 <-> OH_g + *_100',   #22\n\
'OCCH3_{fac} + 3ele_g + 3H2O_g <-> CH3CH2OH_g + *_{fac} + 3OH_g',\n\
#'OCCH3_{fac} + 3ele_g + 2H2O_g <-> C2H4_g + *_{fac} + 3OH_g',\n\
'H2O_g + ele_g + H*_{fac} -> H2_g + *_{fac} + OH_g; beta=0.5',     # Heyrovsky 2\n\
'H*_{fac} + H*_{fac} -> H2_g + 2*_{fac}',     # Tafel 3\n\
#'OH*_{fac} + H_g + ele_g <-> H2O_g + *_{fac}',	#4\n\
#'CH3_{fac} + H2O_g + ele_g <-> CH4_g + *_{fac} + OH_g',	#11\n\
#'CH3_{fac} + H_{fac} <-> CH4_g + 2*_{fac}',	#23\n\
'CO_g + *_{fac} <-> CO_{fac}',\n\
#'C_{fac} + CO_g <-> OCC_{fac}',   #21\n\
]\n".format(fac=facet)
    catmap_string='rxn_expressions=[\n'+catmap_string
    return catmap_string
    #out.write(catmap_string)
#    out.writelines(']')
    #out.close()
    #print(catmap_string)

def collect_field_data(alldata,facet):
    thisfacetdata=field_data[facet].copy()
    #print(thisfacetdata.keys())
    for ads in thisfacetdata.keys():
     if ads in ads_and_electron.keys():
        #print('-'*15)
        charges=[]
        for charge in thisfacetdata[ads].keys():
            if isinstance(charge,float):
                if 'Erel' in thisfacetdata[ads][charge].keys():
                    charges.append(charge)
                else:
                    print(ads+' %s doesnt work'%charge)

        lowest_charge,highest_charge=min(charges),max(charges)
        #XXX: Should be a better fit
        dedq = (thisfacetdata[ads][lowest_charge]['Erel'] - thisfacetdata[ads][highest_charge]['Erel'])/\
                (thisfacetdata[ads][lowest_charge]['qA'] - thisfacetdata[ads][highest_charge]['qA'])

        if 'dedq' not in alldata[ads].keys():
                alldata[ads]['dedq']={}
        if facet not in alldata[ads]['dedq'].keys():
                alldata[ads]['dedq'][facet]=0
        alldata[ads]['dedq'][facet] = dedq
    else:
        print(ads+' seems to be in the field data but not in catmap input')

    for ads in ads_and_electron.keys():
        if ads not in thisfacetdata.keys():
            print(ads+' doesnt seem to have field response data')

def add_base_species_definitions(cat_string,facet):
    cat_string+="\nspecies_definitions['%s'] = {'site_names': ['%s'], 'total':1.0}\n"%(facet,facet)
# site interaction scaling:
    cat_string+="species_definitions['%s']['interaction_response_parameters'] = {'cutoff':0.25,'smoothing':0.05}\n"%facet
    cat_string+="species_definitions['dl'] = {'site_names': ['dl'], 'total':1.0}\n"

    cat_string+="species_definitions['CO_{fac}'] = \
            {{\n\
		'self_interaction_parameter':[None]*numBeforePt+[eCO]+[None]*numAfterPt,\n\
		'cross_interaction_parameters':{{\n\
		'H_{fac}': [None]*numBeforePt+[eHCO]+[None]*numAfterPt,\n\
########	'H2O-ele_t': [None]*numBeforePt+[eHCO*1.08]+[None]*numAfterPt,\n\
########	'H-H2O-ele_t': [None]*numBeforePt+[eHCO*0.7]+[None]*numAfterPt,\n\
########	'H-H_t': [None]*numBeforePt+[eHCO*1.59]+[None]*numAfterPt,\n\
########	'H2O-CO-ele_t': [None]*numBeforePt+[eCO*0.79]+[None]*numAfterPt,\n\
########	'CO-H2O-ele_t': [None]*numBeforePt+[eCO*3.1889]+[None]*numAfterPt,\n\
########	'COH-H2O-ele_t': [None]*numBeforePt+[eCO*0.56]+[None]*numAfterPt,\n\
########	'CH-OH-ele_t': [None]*numBeforePt+[eCO*1.89]+[None]*numAfterPt,\n\
########	'H-CO_t': [None]*numBeforePt+[eCO*1]+[None]*numAfterPt,\n\
########	'OC-CO_t': [None]*numBeforePt+[1.*eCO]+[None]*numAfterPt,\n\
########	'OCCO_t': [None]*numBeforePt+[1.*eCO]+[None]*numAfterPt,\n\
########	'OCCO-OH2-ele_t': [None]*numBeforePt+[1.0238*eCO]+[None]*numAfterPt,\n\
########	'OC-CHO_t': [None]*numBeforePt+[1.194*eCO]+[None]*numAfterPt,\n\
########	'OCCHO_t': [None]*numBeforePt+[1.4839*eCO]+[None]*numAfterPt,\n\
########	'OCCO-H2O-ele_t': [None]*numBeforePt+[1.*eCO]+[None]*numAfterPt,\n\
########	'OCCOH_t': [None]*numBeforePt+[0.44*eCO]+[None]*numAfterPt,\n\
########	'OCC-OH-ele_t': [None]*numBeforePt+[1.*eCO]+[None]*numAfterPt,\n\
########	'OCC_t': [None]*numBeforePt+[1*eCO]+[None]*numAfterPt,\n\
########	'OCC-H_t': [None]*numBeforePt+[1*eCO]+[None]*numAfterPt,\n\
########	'OCC-H2O-ele_t': [None]*numBeforePt+[1*eCO]+[None]*numAfterPt,\n\
########	'OCCH_t': [None]*numBeforePt+[1*eCO]+[None]*numAfterPt,\n\
########	'CO-_t': [None]*numBeforePt+[eCO]+[None]*numAfterPt,\n\
########	'OC-CHOH_t': [None]*numBeforePt+[1.*eCO]+[None]*numAfterPt,\n\
########	'OCCHOH_t': [None]*numBeforePt+[1.2745*eCO]+[None]*numAfterPt,\n\
########	'CHO-OH2-ele_t': [None]*numBeforePt+[0.77*eCO]+[None]*numAfterPt,\n\
########	'CH2O_t': [None]*numBeforePt+[0.3*eCO]+[None]*numAfterPt,\n\
########	'CH2O-H2O-ele_t': [None]*numBeforePt+[0.3*eCO]+[None]*numAfterPt,\n\
########	'CH2OH_t': [None]*numBeforePt+[0.68*eCO]+[None]*numAfterPt,\n\
########	'CH2-OH-ele_t': [None]*numBeforePt+[0.68*eCO]+[None]*numAfterPt,\n\
}} }}".format(fac=facet)

    cat_string+="\n\
species_definitions['OH_{fac}'] = {{\n\
                   'self_interaction_parameter':[None]*numBeforePt+[1.0330]+[None]*numAfterPt,\n\
                   }}\n\n\
species_definitions['H_{fac}'] = {{\n\
                   'cross_interaction_parameters':{{\n\
		    'CO_{fac}': [None]*numBeforePt+[eHCO]+[None]*numAfterPt,\n\
#                   'H2O-ele_t': [None]*numBeforePt+[eH*1.08]+[None]*numAfterPt,\n\
#   #		'H-H2O-ele_t': [None]*numBeforePt+[eH*0.7]+[None]*numAfterPt,\n\
#                   'H-H_t': [None]*numBeforePt+[eH*1.59]+[None]*numAfterPt,\n\
#                   'H2O-CO-ele_t': [None]*numBeforePt+[eHCO*0.79]+[None]*numAfterPt,\n\
#                   'CO-H2O-ele_t': [None]*numBeforePt+[eHCO*3.1889]+[None]*numAfterPt,\n\
#                   'COH-H2O-ele_t': [None]*numBeforePt+[eHCO*0.56]+[None]*numAfterPt,\n\
#                   'H-CO_t': [None]*numBeforePt+[eHCO*1]+[None]*numAfterPt,\n\
                            }}\n\
                   }}\n\
\n\
#   species_definitions['COH_{fac}'] = {{\n\
#                   'self_interaction_parameter':[None]*numBeforePt+[1.75*eCO]+[None]*numAfterPt,\n\
#                   }}\n\
".format(fac=facet)
    cat_string+="\n\
species_definitions['OCCO_{fac}'] = {{\n\
                   'self_interaction_parameter':[None]*numBeforePt+[1*eCO]+[None]*numAfterPt,\n\
                   'cross_interaction_parameters':{{\n\
                   'CO_{fac}': [None]*numBeforePt+[1*eCO]+[None]*numAfterPt,\n\
                   'H_{fac}': [None]*numBeforePt+[1*eHCO]+[None]*numAfterPt,\n\
#                   'H2O-ele_t': [None]*numBeforePt+[1*eHCO*1.08]+[None]*numAfterPt,\n\
#   #		'H-H2O-ele_t': [None]*numBeforePt+[1*eHCO*0.39]+[None]*numAfterPt,\n\
#                   'H-H_t': [None]*numBeforePt+[1*eHCO*1.59]+[None]*numAfterPt,\n\
#                   'H2O-CO-ele_t': [None]*numBeforePt+[1*eCO*0.79]+[None]*numAfterPt,\n\
#                   'CO-H2O-ele_t': [None]*numBeforePt+[1*eCO*3.1889]+[None]*numAfterPt,\n\
#                   'COH-H2O-ele_t': [None]*numBeforePt+[1*eCO*0.56]+[None]*numAfterPt,\n\
#                   'H-CO_t': [None]*numBeforePt+[1*eCO*1]+[None]*numAfterPt,\n\
                           }},\n\
                   'n_sites':2\n\
                   }}\n\
".format(fac=facet)
    return cat_string

def automatize_this_in_the_future():
     for elem in ['C','O','H']:
          elem_ads[elem] = ads.replace('H2','HH').replace('H3','HHH').count(elem)#unique(ads,return_counts=True))
          elem_toads[elem] = toads.replace('H2','HH').replace('H3','HHH').count(elem)#unique(ads,return_counts=True))

     #XXX: More than elementary steps are not implemented atm
     if ads_and_electron[toads_in]['nHe'] - ads_and_electron[ads_in]['nHe'] > 1:
                print(ads_in+' to '+toads_in+' doesnt seem to be an elementary step')
                #continue


     if elem_ads['C'] == elem_toads['C']:
         if elem_toads['O'] < elem_ads['O']:
             if elem_toads['H'] < elem_ads['H']:
                 #OH desorbing
                 catmap_line += OH_des
             else:
                 print('Weird! From %s to %s an O left without prior protonation'%(ads,toads))
                 #continue
         elif elem_toads['O'] == elem_ads['O']:
             if elem_toads['H'] > elem_ads['H']:
                 #Protonation
                 catmap_line += O_H_ads#'ele_g + H2O_g <-> %s_%s + OH_g'%(toads,facet)
             else:
                 print('Weird! From %s to %s an oxidation is  happening which is  not  implemented yet')
                 #continue

         if ads_and_electron[ads_in]['n_sites'] != ads_and_electron[toads_in]['n_sites']:
#                 print('From %s to %s the dentation changed'%(ads_in,toads_in))
                 ns_ads = ads_and_electron[ads_in]['n_sites']
                 ns_toads = ads_and_electron[toads_in]['n_sites']

                 if ns_ads  > ns_toads:
                     catmap_line+=' + *_%s'%site[toads_in]
                 else:
                     catmap_line='*_%s + '%site[ads_in]+catmap_line

     else:
                ads2 = toads.replace(ads,'',1)
                if ads2 not in ads_and_electron.keys():
                    if ads2[::-1] not in ads_and_electron.keys():
                        print('Creating dummy at %s to %s, since the second reactant can not be created'%(ads,toads))
                        #catmap_line+='%s_%s <-> %s_%s'%('XXX',facet,toads,facet)
                        ads_str='XXX_%s'%site[ads_in]
                    else:
                        ads_str=ads2[::-1]+'_'+site[ads_in]
                        #catmap_line+='%s_%s <-> %s_%s'%(ads2[::-1],facet,toads,facet)
                else:
                    if ads2 == toads:
                        ads_str=''
                        #if barriers is None:
                        #    catmap_line+='<-> %s_%s'%(toads,site[toads_in])
                        #else:
                        #    catmap_line+='<-> ^%1.4feV_%s <-> %s_%s'%(barriers['C-C'],site[ads_in],toads,site[toads_in])
                    else:
                        ads_str=ads2+'_'+site[ads_in]
                        #catmap_line+='%s_%s <-> %s_%s'%(ads2,facet,toads,facet)

                if barriers is None:
                        catmap_line+='%s <-> %s_%s'%(ads_str,toads,site[toads_in])
                else:
                    if 'E_ddag_%s'%facet in ads_and_electron[ads_in].keys():
                        for pH in ['base','acid']:
                        #print(ads,toads,ads_and_electron[ads_in]['E_ddag_%s'%facet])
                        #print(ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet])
                            TSname=get_unique_TS_name(ads,toads,pH,written_barriers)
                        #print(ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet])
                        catmap_line+='%s <-> %s_%s <-> %s_%s; beta=%s'%(
                                ads_str,TSname,site[ads_in],toads,site[toads_in],
                            ads_and_electron[ads_in]['E_ddag_vs_pot_%s'%facet][list(ads_and_electron[ads_in]['E_ddag_%s'%facet].keys())[0]][0])
                        written_barriers.append(TSname)

                    else:
                        catmap_line+='%s <-> ^%1.4feV_%s <-> %s_%s'%(ads_str,barriers['C-C'],site[ads_in],toads,site[toads_in])

                if ads_and_electron[ads_in]['n_sites'] == ads_and_electron[toads_in]['n_sites']:
                    catmap_line+=' + *_%s'%facet
