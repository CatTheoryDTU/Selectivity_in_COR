import matplotlib
import pickle
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import os,sys
import numpy as np
from general_tools import get_reference_energies,lin_fun,quad_fun,get_reference_vibrational_contribution
from ase import units
from sjm_analyse_tools import get_intcolors

def plot_pseudo_endstate_difference_to_real_endstate(dat,VSHE=4.4,pH=13,entype='E',facet='100'):
    outdat=[]
    names=[]
    counter=0
    for ads in dat:
        baren=entype+'_ddag_vs_pot_'+facet
        if baren not in dat[ads]: continue
        for toads in dat[ads][baren]:
            #for pHnames in ['acid','base']: #chemical does not have a pseudo endstate
            for pHnames in ['base']: #chemical does not have a pseudo endstate
                if pHnames not in dat[ads][baren][toads]: continue

                IS = dat[ads][entype+'_vs_pot_100']
                TS = dat[ads][entype+'_ddag_vs_pot_100'][toads][pHnames]
                if pHnames == 'base':
                #    if 'base' not in dat[ads][entype+'_ddag_FS_vs_pot_100'][toads]: continue
                    pES = dat[ads][entype+'_ddag_FS_vs_pot_100'][toads]['base']
                else:
                    if 'acid' not in dat[ads][entype+'_ddag_IS_vs_pot_100'][toads]: continue
                    pES = dat[ads][entype+'_ddag_IS_vs_pot_100'][toads]['acid']
                print(toads)
#    pIS = dat[ISname][entype+'_ddag_IS_vs_pot_100'][FSname]['base']

                if toads in ['CCO','HCCO','CCH2','CCHO']:
                    FS = dat['md-'+toads][entype+'_vs_pot_100']
                elif toads in ['HCCHO','H2CCH2O']:
                    FS = dat['bdo-'+toads][entype+'_vs_pot_100']
                else:
                    try:
                        FS = dat[toads][entype+'_vs_pot_100']
                    except:
                        print(toads)
                        continue

                FS[1]=FS[1]-(4.4-0.059*pH)
                FS[0]+=1
                if pHnames == 'acid':
                    outdat.append([counter,(pES[0]-IS[0])*VSHE+(pES[1]-IS[1])])
                else:
                    outdat.append([counter,(FS[0]-pES[0])*VSHE+(FS[1]-pES[1])])
                    #print(ads,toads,pHnames)
                names.append(get_intcolors(ads,toads,return_name=True)+pHnames[0])
                counter+=1
    outdat=np.array(outdat)
    plt.plot(outdat[:,0],outdat[:,1],'o')

    plt.xticks(range(len(names)),names,rotation=90)
#plt.savefig('Band_with_pFS.pdf')
    plt.tight_layout()
    plt.show()


