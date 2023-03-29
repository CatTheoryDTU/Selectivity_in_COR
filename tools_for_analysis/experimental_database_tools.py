#from plot_paper_data import load_pkl_data, _compute_partial_jecsa
import matplotlib.patheffects as path_effects
from matplotlib import pyplot as plt
#from scipy.optimize import curve_fit
from ase import units
import numpy as np
import string

#!/usr/bin/env python

import sys, os
import numpy as np
from scipy.optimize import curve_fit, leastsq
from copy import deepcopy
import pickle

#import matplotlib.pyplot as plt
#from rtools.helpers.matplotlibhelpers import tumcolors as tumcs
#import rtools.helpers.matplotlibhelpers as matplotlibhelpers
from scipy.constants import golden_ratio, inch
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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

def _set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.5/2,\
                    lrbt=[0.135,0.80,0.25,0.95],fsize=9.0,font='helvetica'):
    # set plot geometry
    rcParams['figure.figsize'] = (width, height) # x,y
    rcParams['font.size'] = fsize
    rcParams['figure.subplot.left'] = lrbt[0]  # the left side of the subplots of the figure
    rcParams['figure.subplot.right'] = lrbt[1] #0.965 # the right side of the subplots of the figure
    rcParams['figure.subplot.bottom'] = lrbt[2] # the bottom of the subplots of the figure
    rcParams['figure.subplot.top'] = lrbt[3] # the bottom of the subplots of the figure

    rcParams['xtick.top'] = True
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.right'] = True
    rcParams['ytick.direction'] = 'in'

    rcParams['legend.fancybox'] = False
    #rcParams['legend.framealpha'] = 1.0
    rcParams['legend.edgecolor'] = 'k'
    if font != 'None':
        matplotlibhelpers.set_latex(rcParams,font=font) #poster

def _compute_partial_jecsa(data, adsorbate, maxpH=100.,minpH=0):
    data = deepcopy(data)
    out = {}
    for k in data:
        if data[k]['roughness'] != 'None' and adsorbate in data[k]:
            # transform total current into j-ecsa
            try:
                pj = _partial_jecsa(data[k], adsorbate)
            except:
                continue
            else:
                ind = np.where(pj != 0.0)[0]
                pj = np.absolute(pj[ind,:]) # correct for current definition

                urhe = data[k]['V_RHE'][ind,:]
                ushe = deepcopy(urhe)
                ushe[:,0] -= 0.059 * float(data[k]['pH'])
            # multiply with FE for chosen product

        #    pj = np.array([data[k][adsorbate][i,:]*j_ecsa[i] for i in range(len(j_ecsa))])
        #    ind = np.where(pj != 0.0)[0]
        #    pj = np.absolute(pj[ind,:]) # correct for current definition
        #    urhe = data[k]['V_RHE'][ind,:]
        #    ushe = deepcopy(urhe)
        #    ushe[:,0] -= 0.059 * float(data[k]['pH'])
            if minpH <= float(data[k]['pH']) <= maxpH:
                out.update({k:{'U_RHE':urhe, 'U_SHE':ushe, 'j_partial':pj, 'cell':data[k]['cell'], \
                    'mode':data[k]['mode'], 'pH':data[k]['pH'], 'roughness':data[k]['roughness'],
                    'folder':data[k]['folder']}})
        else:
            print("%s does not contain roughness or %s"%(k, adsorbate))
    return(out)

def _partial_jecsa(dat, adsorbate):
    j_ecsa = dat['I(mA/cm2)'] / float(dat['roughness'])
    # multiply with FE for chosen product
    pj = np.array([dat[adsorbate][i,:]*j_ecsa[i] for i in range(len(j_ecsa))])
    return(pj)

# global keywords
c2a = ['Ethylene', 'Acetate', 'Ethanol', 'Acetaldehyde', 'n-propanol',\
        'Allyl-Alcohol','Acetone','Ethane','Propionaldehyde','Allyl-alcohol',\
        'Hydroxyacetone','Ethylene-Glycol','Glycolaldehyde','Ethylene-glycol']
c1a = ['Methane', 'Methanol', 'Formate']


def _add_c1_c2(data):
    for k in data:
        # check for new adsorbates
        for e in data[k]:
            if e not in c2a+c1a+['Hydrogen','V_RHE','mode','pH','catalyst','roughness','cell','doi','I(mA/cm2)','CO']:
                print(e)
        # sum FE of C1
        c1 = np.zeros((data[k]['V_RHE'].shape))
        for e in c1a:
            if e in data[k]:
                c1[:,0] += data[k][e][:,0]
                c1[:,1] += data[k][e][:,1]**2.0
        c1[:,1] = np.sqrt(c1[:,1])
        data[k].update({"C1":c1})
        # sum FE of C2s
        c2 = np.zeros((data[k]['V_RHE'].shape))
        for e in c2a:
            if e in data[k]:
                c2[:,0] += data[k][e][:,0]
                c2[:,1] += data[k][e][:,1]**2.0
        c2[:,1] = np.sqrt(c2[:,1])
        data[k].update({"C2+":c2})
    return(data)

def _load_pickle_file(filename,py23=False):
    pickle_file = open(filename,'rb')
    if py23:
        data = pickle.load(pickle_file, encoding='latin1')
    else:
        data = pickle.load(pickle_file)
    pickle_file.close()
    return(data)

def _write_pickle_file(filename,data,py2=False):
    ptcl = 0
    if py2:
        filename = filename.split('.')[0]+'_py2.'+filename.split('.')[1]
        ptcl = 2
    output = open(filename, 'wb')
    pickle.dump(data, output, protocol=ptcl)
    output.close()

def plot_partial_current_densities(filename, data, pot, clr='data'):
    _set_plotting_env(width=3.37*1.3,height=3.37,\
                   lrbt=[0.15,0.75,0.13,0.98],fsize=7.0,font='None')

   #vcolors = [tumcs['tumorange'],tumcs['tumgreen'],tumcs['tumblue'],\
   #    tumcs['tumlightblue'],tumcs['acc_yellow'],tumcs['diag_violet']]
    if clr == 'data':
        ks = list(data.keys())
    elif clr == 'pH':
        ks = (np.unique([float(data[k]['pH']) for k in data])).tolist()
    colors = plt.cm.jet(np.linspace(0,1,len(ks)))
    kclr = {ks[i]:colors[i] for i in range(len(ks))}
    lss = {'COR':'-', 'CO2R':'--'}
    mks = {'H-cell':'o', 'GDE':'x'}


    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    for k in data:
        x1 = data[k][pot][:,0]; xerr1 = data[k][pot][:,1]
        #x2 = data[k]['U_RHE'][:,0]; xerr2 = data[k]['U_RHE'][:,1]
        if np.all(xerr1 == 0.0):
            xerr1 = None
       #if np.all(xerr2 == 0.0):
       #    xerr2 = None
        kc = k
        if clr == 'pH':
            kc = float(data[k]['pH'])

        y = data[k]['j_partial'][:,0]; yerr = data[k]['j_partial'][:,1]
        if np.all(yerr == 0.0):
            yerr = None
        #if True:
        if data[k]['cell'] == 'GDE' and data[k]['mode'] == 'COR':
            ax1.errorbar(x1, y, xerr=xerr1, yerr=yerr, color=kclr[kc], ls=lss[data[k]['mode']], marker=mks[data[k]['cell']], markersize=3)
            ax2.errorbar(x1, y, xerr=xerr1, yerr=yerr, color=kclr[kc], ls=lss[data[k]['mode']], marker=mks[data[k]['cell']], markersize=3)

    # legend
    for l in lss:
        ax1.plot(np.nan, np.nan, color='k', ls=lss[l], label=r'%s'%l)
    for m in mks:
        ax1.plot(np.nan, np.nan, color='k', ls='None', marker=mks[m], label=r'%s'%m)
    for k in ks:
        if data[k]['cell'] == 'GDE' and data[k]['mode'] == 'COR':
            ax1.plot(np.nan, np.nan, color=kclr[k], label=r'%s'%k)
    ax1.legend(loc=1,prop={'size':4},bbox_to_anchor=(1.45, 1.))
    ax2.set_zorder(-1)
    #ax2.set_ylim(1e-5,10)
    ax2.set_ylim(1e-5,100)

    # axis labels
    ax1.set_xticklabels([])
    ax2.set_yscale('log')
    ax2.set_xlabel('U vs. %s (V)'%pot.split('_')[1])
    ax2.set_ylabel('$j_{\mathrm{ECSA}}$ (mA/cm$^2$)')
    ax2.yaxis.set_label_coords(-0.15, 1.05)
    plt.subplots_adjust(hspace=0.05)
    #plt.show()
    #matplotlibhelpers.write(filename+'_'+pot,folder='output',\
#        write_info=False,write_png=False,write_pdf=True,write_eps=False)
    writefig(filename+'_'+pot,folder='output')


def writefig(filename, folder='output',  write_eps=False):
    """
      wrapper for creating figures
      Parameters
      ----------
      filename : string
        name of the produced figure (without extention)
      folder : string
        subfolder in which to write the figure (default = output)
      write_eps : bool
        whether to create an eps figure (+ the usual pdf)
    """
    # folder for output
    if not os.path.isdir(folder):
        os.makedirs(folder)
    fileloc = os.path.join(folder, filename)
    print("writing {}".format(fileloc+'.pdf'))
    plt.savefig(fileloc+'.pdf')
    if write_eps:
        print("writing {}".format(fileloc+'.eps'))
        plt.savefig(fileloc+'.eps')
    plt.close(plt.gcf())

#def create_input_pkl(folder):
#    scripts = [ script for script in os.listdir(folder)
#            if script[-3:] == '.py']
#    basedir=os.getcwd()
#    if len(scripts) > 1:
#        print(CRED+"Careful! There's more than one python script in %s"%folder+CEND)
#    #os.chdir(basedir+'/'+folder)
#    os.system('python3 %s'%scripts[0])
    #os.chdir(basedir)

def load_pkl_data(folders,data=None):
    if data is None:
        data={}
    # collect all data
    for f in folders:
        if not os.path.exists(f) or not os.path.isdir(f):
            print(CRED+'Couldnt find the folder %s, check the spelling!'%f+CEND)
            continue

        pklfile = "%s/%s.pkl"%(f,f.split('/')[-1])
        if pklfile.split('/')[-1] not in os.listdir(f):
            print('\033[92mCouldnt find input pkl in %s, creating it\033[0m'%f)
            prep_exp_data(f)

        dat = _load_pickle_file(pklfile)

        tag = f.split('_')[0][:5]
        for k in dat:
            #print(k)

            if k[:5] != 'Cu-oh': # faulty potential
                data.update({tag+'_'+k:dat[k]})
                data[tag+'_'+k]['folder']=f
    return data


def _read_header(txtfile):
    with open(txtfile, 'r') as infile:
        hline = infile.readlines()[0][1:]
    return(hline)

def _read_meta_data(mfile):
    with open(mfile, 'r') as infile:
        lines = infile.readlines()
    dat = {}
    for l in lines:
        e = [s.strip() for s in l.split('=')]
        dat.update({e[0]:e[1]})
    return(dat)

def prep_exp_data(study):
    # read main data
    # no std std deviation
    full_data2 = {}
    nfiles,mfiles=[],[]
    #study = _check_execution_dir(f)
    #print(study)
    for files in os.listdir(study):
        if files.split('_')[0] == 'raw':
            if files.split('_')[1][:4] == 'data':
                nfiles.append(study+'/'+files)
        elif files.split('_')[0] == 'meta':
            if files.split('_')[1][:4] == 'data':
                mfiles.append(study+'/'+files)
    if len(nfiles) != len(mfiles):
        print('The number of experimental and meta data files in %s do not\
               match!  Check!'%(os.getcwd().split('/')[-1]))


    for nfile in sorted(nfiles):
        #nfile = 'raw_data_%i.txt'%i
        keys = _read_header(nfile).split()
        val = np.loadtxt(nfile)
        if len(val.shape) == 1:
            val=np.array([val])

        if 'I(mA/cm2)' not in keys:
            #print(study)
            if 'FE_tot' in keys:
                # scaled by FE tot and FE of products
               jsum = val[:,2:].sum(axis=1)
               jtot = jsum/(val[:,1]/100.)
               for j in range(len(jtot)):
                   val[j,2:] /= jtot[j]
               data2 = {keys[k]:np.zeros((val.shape[0], 2)) for k in range(len(keys)) if keys[k] != 'FE_tot'}
               for k in data2:
                  data2[k][:,0] = val[:,keys.index(k)]
               v = np.zeros((len(jtot),2)); v[:,0] = jtot
               data2.update({'I(mA/cm2)':v})
            else:
               continue
        else:
        #if  1:
            val[:,2:] /= 100. # correct for percentage

            # prep data
            data2 = {keys[i]:np.zeros((val.shape[0], 2)) for i in range(len(keys))}
            for k in data2:
                data2[k][:,0] = val[:,keys.index(k)]

        mfile = '/'.join([nfile.split('/')[0],'meta'+nfile.split('/')[1].lstrip('raw')])
        #print(nfile,mfile,os.listdir())
        #print(data2)
        if mfile.split('/')[-1] in os.listdir(study):
            mdat = _read_meta_data(mfile)
            data2.update(mdat)
            data2['folder']=study

            # sort into one file
            full_data2.update({mdat['catalyst']+'_pH'+str(mdat['pH']):data2})

    #print(full_data2)
    #oo
    # save data
    #name = os.getcwd().split('/')[-1]
    _write_pickle_file(study+"/%s.pkl"%study, full_data2)


def main():
    if 1:
        fig,ax=plt.subplots(1,1,figsize=(10,7))
        #data=read_data(folders,emphasized_folders)
        #studies = get_studies(products,data,100,0,[],False)
        studies,folderlabels,all_folders,fits=perform_Tafel_analysis(folders,
#                ['Ethylene+Ethanol+Acetate+n-propanol'],
                ['+'.join(products)],
                **general_kw,
                verbose=True,
                ax=ax)
        plt.legend(bbox_to_anchor=(1.,1),fontsize=8.0,loc='upper left')
        plt.tight_layout()
        plt.savefig('output/Tafel_analysis_all_products.pdf')
        #plt.show()
        plt.close()
        #dd


    bakk_kw=general_kw['mode']
    general_kw['mode'] = 'COR'
    d1,d2,d3,fitsCO = perform_Tafel_analysis(folders,
            products,
            **general_kw,
            studies=studies,folderlabels=folderlabels)
    general_kw['mode'] = 'CO2R'
    d1,d2,d3,fitsCO2 = perform_Tafel_analysis(folders,
            products,
            **general_kw,
            studies=studies,folderlabels=folderlabels)
    plt.close()
    general_kw['mode']=bakk_kw

    offset=0.08
    axw=0.44

    fig=plt.figure(figsize=(12,12))
    axSCO2= fig.add_axes([offset,offset,axw,0.4])
    axSCO = fig.add_axes([offset+axw,offset,axw,0.4])
    S = plot_selectivity(fitsCO,colors=['b','peru','r','g','k'],products=['Ethylene','Ethanol','Acetate','n-propanol','Ethanol+Acetate+n-propanol'],plot_FE=False,ax=axSCO,mode='COR')
    S = plot_selectivity(fitsCO2,colors=['b','peru','r','g','k'],products=['Ethylene','Ethanol','Acetate','n-propanol','Ethanol+Acetate+n-propanol'],plot_FE=False,ax=axSCO2,mode='CO2R')
    #S = plot_selectivity(fitsCO,colors=['b','peru','r','g','k'],products=['Ethylene','Ethanol','Acetate','Ethanol+Acetate'],plot_FE=False,ax=axSCO,mode='COR')
    #S = plot_selectivity(fitsCO2,colors=['b','peru','r','g','k'],products=['Ethylene','Ethanol','Acetate','Ethanol+Acetate'],plot_FE=False,ax=axSCO2,mode='CO2R')
    axSCO.set_ylabel('')
    axSCO.set_yticks([])
    axSCO2.annotate('CO$_2$R',(-1.4,0.97), fontsize=24,fontweight='bold',color='k',va='top').draggable()
    axSCO.annotate('COR',(-1.4,0.97), fontsize=24,fontweight='bold',color='k',va='top').draggable()
    axSCO2.annotate('(e)',(-1.68,0.99),va='top',ha='left',fontsize=20)
    axSCO.annotate('(f)',(-1.68,0.99),va='top',ha='left',fontsize=20)

    axes=[]
    axh=0.2
    lowb=0.54
    highb=0.76
    axes.append(fig.add_axes([offset,lowb+axh,axw,axh]))
    axes.append(fig.add_axes([offset,lowb,axw,axh]))
    axes.append(fig.add_axes([axw+offset,lowb+axh,axw,axh]))
    axes.append(fig.add_axes([axw+offset,lowb,axw,axh]))

    for iprod,product in enumerate(products):
          perform_Tafel_analysis(folders,
                product,
                **general_kw,
                ax=axes[iprod],
                studies=studies,folderlabels=folderlabels)
          tag='abcd'[iprod]
          axes[iprod].annotate(f'({tag})',(-1.68,-4.2),va='center',ha='left',fontsize=20)


    plt.savefig('output/Fig4.pdf')

    plt.show()



def plot_selectivity(fits,plot_FE=True,products=['Ethylene','Ethanol+Acetate+n-propanol'],ax=None,pots=None,colors=None,outdir='output',add_prod_label=True,
        linestyle='-',mode=None):

    if colors is None:
        colors=plt.cm.RdYlGn(np.linspace(0, 1, len(products)))
    if ax is None: fig,ax = plt.subplots()
    if pots is None: pots=np.arange(-1.7,-0.94,0.05)

    elec={'Ethylene':8,'Ethanol':8,'Acetate':4,'Acetaldehyde':6,'n-propanol':12}
    if mode=='COR':
        elec={'Ethylene':8,'Ethanol':8,'Acetate':4,'Acetaldehyde':6,'n-propanol':12}
    elif mode == 'CO2R':
        elec={'Ethylene':12,'Ethanol':12,'Acetate':8,'Acetaldehyde':10,'n-propanol':18}
    else:
        print('Using the number of electrons from COR, but be careful it is not correct for CO2R')

    totfit=np.array([0.,0.])
    #print(fits)
    for prod in fits:
        totfit+=fits[prod]
    totatpots=[]
    fitsatpots={}

    #Get the total current (or rate if plot_FE=False)
    totatpots = np.zeros((len(pots)))#10**(totfit[0]*pots+totfit[1])
    for prod in fits:
        fitsatpots[prod] = 10**(fits[prod][0]*pots+fits[prod][1])
        if not plot_FE :
            fitsatpots[prod] /= elec[prod]
        totatpots+=fitsatpots[prod]

    #Calculate FE or selectivities (if plot_FE=False)
    selectivities={}
    for prod in fits:
        selectivities[prod]=[]
        for ipot in range(len(pots)):
            selectivities[prod].append([pots[ipot],fitsatpots[prod][ipot]/totatpots[ipot]])
        selectivities[prod]=np.array(selectivities[prod])
    #    plt.plot(selectivities[prod][:,0],selectivities[prod][:,1],label=prod)


    #Plot selectivities based "products"
    for iprod,prod in enumerate(products):
        prods=prod.split('+')
        sel=None

        for singleprod in prods:
            if sel is None:
                sel=selectivities[singleprod][:,1]
            else:
                sel+=selectivities[singleprod][:,1]

    #        if xvals is None:
    #            xvals=selectivities[singleprod][:,0]
        ax.plot(pots,sel,color=colors[iprod],linewidth=2,linestyle=linestyle)
        if add_prod_label:
            if len(prods) == 1:
                ax.annotate(prod,(pots[-1],sel[-1]),ha='right',color=colors[iprod],fontsize=24).draggable()
            else:
                ax.annotate('Oxygenates',(pots[-1],sel[-1]),ha='right',color=colors[iprod],fontsize=24).draggable()
        #    prod='Ethylene'
#    ax.plot(selectivities[prod][:,0],selectivities[prod][:,1],'k')
#    ac_eoth_sel=[]


#    for ipot in range(len(pots)):
#        ac_eoth_sel.append(selectivities['Ethanol'][ipot][1]+selectivities['Acetate'][ipot][1]+selectivities['n-propanol'][ipot][1])
#    ac_eoth_sel=np.array(ac_eoth_sel)

#    plt.plot(selectivities[prod][:,0],ac_eoth_sel,'r')
    ax.set_ylim([0,1])
    ax.set_xlim([-1.69,-0.95])
    ax.set_xlabel('U$_{\mathrm{SHE}}$ / V')
    if plot_FE:
        ax.set_ylabel('Faradaic efficiency / %')
    else:
        ax.set_ylabel('Selectivity')

    if plot_FE:
        plt.savefig(outdir+'/Faradaic_efficiencies.pdf')
    else:
        plt.savefig(outdir+'/Selectivities.pdf')
    #plt.legend()
#    plt.show()
    return selectivities


def compare_products_in_study(folders,products,maxpH=100,potrange=None,
        potscale='SHE',reference=None,outdir='output',max_standard_deviation=0.8):
    c2ads=['Ethylene','Ethanol','Acetate','n-propanol']
    comolads=['Ethylene','Ethanol','Acetate','n-propanol','Methane']
    if potrange is None:
        potrange=np.arange(-2,-1,0.1)

    if 'CO_mol' in products:
        COcons=True
    else:
        COcons=False

    fits = plot_totalcurrents(folders,[],products,
        maxpH=100,potrange=potrange,potscale=potscale,max_standard_deviation=max_standard_deviation,
        COconsumption=COcons,outdir=outdir,plot_results=False)

    if reference is None:
        reference=list(fits.keys())[-1]

    refcurrents={}
    for totads in products:
      refcurrents[totads]=0
      for ads in totads.split('+'):
          if ads == 'C2':
             for c2s in c2ads:
                 refcurrents[totads]+=fits[reference][c2s]
          elif ads == 'CO_mol':
             for c2s in comolads:
                 refcurrents[totads]+=fits[reference][c2s]
          else:
              refcurrents[totads]+=fits[reference][ads]

    relcurrents={}
    for istudy, study in enumerate(fits.keys()):
        if study == reference: continue
        relcurrents[study]={}
        for totads in products:
            relcurrents[study][totads]=0
            for ads in totads.split('+'):
                if ads == 'C2':
                    for c2s in c2ads:
                        relcurrents[study][totads]+=fits[study][c2s]
                elif ads == 'CO_mol':
                    for c2s in comolads:
                        relcurrents[study][totads]+=fits[study][c2s]
                else:
                    relcurrents[study][totads]+=fits[study][ads]
            relcurrents[study][totads]/=refcurrents[totads]

    colors = plt.cm.jet(np.linspace(0, 1, len(relcurrents.keys())))
    markers={'Methane':'^',
            'Ethylene': 'o',
            'Ethanol': 'P',
            'n-propanol': '>',
            'Acetate': '<',
            'Hydrogen':'d',
            'C2':'h',
            'CO_mol': '*'
            }
    plt.axhline(y=1,linestyle='--',color='k')
    for istudy, study in enumerate(relcurrents.keys()):
        plt.plot(np.nan,np.nan,'s',color=colors[istudy],label=study)
        for totads in relcurrents[study].keys():
            plt.plot(fits[study]['V_SHE'],relcurrents[study][totads],markers[totads],color=colors[istudy],markeredgecolor='k')
    for totads in relcurrents[study].keys():
        plt.plot(np.nan,np.nan,markers[totads],color='k',label=totads)

    plt.legend(bbox_to_anchor=(1,1))
    plt.ylabel('Relative current vs %s'%reference)
    plt.xlabel('U_\mathrm{%s} / V'%potscale)
    plt.tight_layout()
    outname = 'Reative_current_%s_%s_maxpH%s_U_\mathrm{%s}'% ('-'.join(folders),'-'.join(products),maxpH, potscale)
    plt.savefig(outdir+'/'+outname+'.pdf')
    #plt.show()
    plt.close()


def plot_totalcurrents(folders,emphasized_folders,products,
        maxpH=100,potrange=None,potscale='SHE',max_standard_deviation=0.98,
        COconsumption=True,outdir='output',plot_results=True):

    colors = ['gray','r','b','g','y']#plt.cm.jet(np.linspace(0, 1, len(a_dat.keys())))
    if potrange is None:
        potrange=np.arange(-2,-1,0.1)
    data = read_data(folders,emphasized_folders)

    allproducts=[]
    nelecs={'Methane':6,
            'Ethylene': 8,
            'Ethanol': 8,
            'n-propanol': 8,
            'Acetate': 4,
#            'Hydrogen':2
            }
    totads=[]
    fitted_data={}
    for ads in nelecs.keys():

#        if '+' in ads: continue
        a_dat = _compute_partial_jecsa(data, adsorbate=ads, maxpH=maxpH)
        for study in a_dat.keys():
            if study not in fitted_data.keys():
                fitted_data[study]={'V_%s'%potscale: potrange}

            if not all(par in a_dat[study].keys()
                   for par in ['U_'+potscale, 'j_partial']):
                print(CRED + 'Couldnt find U_%s and j_partial in %s'
                      % (potscale, study) + CEND)
                continue

            # Perform Tafel fit
            coeff = preselect_datapoints(a_dat[study],
                                     max_standard_deviation,
                                    potscale,ads,add_pH_to_rate=add_pH_to_rate)[2]
            #if ads == 'Methane':
                #print(study,a_dat[study]['U_SHE'][:,0])
                #print(study,a_dat[study]['j_partial'][:,0])
            # Fitting current from Tafel analysis so potentials are equal
            # for all products
            if coeff is None: continue
            j_fit=lin_fun(potrange,*coeff)

            #Changing units to mol/(cm^2*s)
            if COconsumption:
                mol_fit=10**j_fit*units.C/(1e3*units.mol*nelecs[ads])
                if ads in ['Ethylene','Ethanol','Acetate']:
                    mol_fit*=2
                elif ads in ['n-propanol']:
                    mol_fit*=3

                fitted_data[study][ads]=mol_fit
            else:
                fitted_data[study][ads]=10**j_fit

    for istudy,study in enumerate(fitted_data.keys()):
        summol=None
        for ads in fitted_data[study].keys():
            if (ads not in nelecs.keys() or
                ads in ['Hydrogen']): continue

            if summol is None:
                summol = fitted_data[study][ads].copy()
            else:
                summol += fitted_data[study][ads].copy()
        fitted_data[study]['total_CO_moles']=summol


        if plot_results:
            plt.plot(fitted_data[study]['V_%s'%potscale],fitted_data[study]['total_CO_moles'],color=colors[istudy],label='%s: %1.1f mV/dec'                             % (study, 1000/coeff[0]))

    if not plot_results:
        return fitted_data

    plt.yscale('log')
    plt.xlabel(r'V$_{\mathrm{SHE}}$ [V]')
    plt.legend(bbox_to_anchor=(1, 1),prop={'size':4})
    if COconsumption:
        plt.ylabel(r'CO comsumption [$\frac{mol}{s cm^2}$]')
        outname=outdir+'/CO_molconsumption_maxpH%s_U_\mathrm{%s}'% (maxpH, potscale)
    else:
        plt.ylabel(r'Total current [$\frac{mA}{cm^2}$]')
        outname=outdir+'/Total_current_maxpH%s_U_\mathrm{%s}'% (maxpH, potscale)
    name='CO_consumption_'
    plt.tight_layout()

    plt.savefig(outname+'.pdf')
    plt.close()
    return fitted_data


def get_studies(products,data,maxpH,minpH,emphasized_folders,color_by_pH):
    studies={}
    for totads in products:
      for ads in totads.split('+'):
        a_dat = _compute_partial_jecsa(data, adsorbate=ads, maxpH=maxpH,minpH=minpH)
        for istudy, study in enumerate(a_dat.keys()):
            if study not in  studies:
                studies[study]={}
    if len(emphasized_folders):
        nemphstudies=0
        for study in data:
            if 'emphasized' in data[study]:
                nemphstudies+=1
        colors = plt.cm.jet(np.linspace(0, 1, nemphstudies))
    elif color_by_pH:
        pHs=np.unique([np.around(float(data[study]['pH']),2) for study in data])
        colors= plt.cm.Blues(np.linspace(0,1,len(pHs)))
    else:
        colors=plt.cm.jet(np.linspace(0, 1, len(studies.keys())))

    iemps=0
    for istudy, study in enumerate(studies.keys()):
        if len(emphasized_folders):
            if 'emphasized' in data[study]:
                studies[study]['color']=colors[iemps]
                iemps+=1
            else:
                studies[study]['color']='lightgray'
        elif color_by_pH:
            if study in data: #len(studies[study]):
                studies[study]['color']=colors[np.where(np.around(float(data[study]['pH']),2) == np.array(pHs))[0][0]]
                #print(studies[study]['color'])
            else:
                #print(study,a_dat)
                studies[study]['color']='r'
        else:
                studies[study]['color']=colors[istudy]

    return studies

def perform_Tafel_analysis(folders, products, emphasized_folders=[], maxpH=100, outdir='output',
                           min_datapoints=2, max_standard_deviation=0.98, minpH=0,xlim=[-1.69,-0.61],
                           mode=None, potscale='SHE',plot_individual_fits=False,individual_Tafel_slopes_in_label=False,
                           fit_CO_and_CO2_separately=False,add_pH_to_rate=False,studies=None,folderlabels=None,ax=None,
                           color_by_pH=False,verbose=False,show_rsquared=False):
    """
    Main function for making the plots
    """
    if ax is None:
        fig,ax=plt.subplots()
    data=read_data(folders,emphasized_folders)
    pot_response_fits={}
    nfolders=0
    if folderlabels == None:
        folderlabels={}
    markers={'CO2R':'d','COR':'o'}
    if studies==None:
#        studies={}
        studies = get_studies(products,data,maxpH,minpH,emphasized_folders,color_by_pH)


    if isinstance(products,str):
        products=[products]

    for totads in products:
      print('Tafel slopes for %s in studies:' % totads)
      tot_all_data, Tafdata, Tafred,all_folders,all_modes = {}, [], [],{},{}
      all_labels=[]

      #plt.plot(np.nan,np.nan,'d',markeredgecolor='k',markerfacecolor='w',label='CO$_2$R')
      #plt.plot(np.nan,np.nan,'o',markeredgecolor='k',markerfacecolor='w',label='COR')
      #plt.plot(np.nan,np.nan,'d',markeredgecolor='w',markerfacecolor='w',label='-'*20)
      COR_Tafred,CO2R_Tafred=[],[]
      for ads in totads.split('+'):
        a_dat = _compute_partial_jecsa(data, adsorbate=ads, maxpH=maxpH,minpH=minpH)
        for istudy, study in enumerate(a_dat.keys()):
            folder=a_dat[study]['folder']
            if folder not in folderlabels:
                    folderlabels[folder]=string.ascii_uppercase[nfolders]
                    nfolders+=1
            if mode and a_dat[study]['mode'] != mode:
                continue
            if study not in tot_all_data.keys():
                tot_all_data[study] = []

            if not all(par in a_dat[study].keys()
                       for par in ['U_'+potscale, 'j_partial']):
                print(CRED + 'Couldnt find U_%s and j_partial in %s'
                      % (potscale, study) + CEND)
                continue
            folder=a_dat[study]['folder']
            all_folders[study]=folder
            all_modes[study]=a_dat[study]['mode']
            this_reddata, this_all_data, coeff = \
                preselect_datapoints(a_dat[study],
                                     max_standard_deviation,
                                     potscale,ads,add_pH_to_rate=add_pH_to_rate)

            alpha=1 if 'emphasized' in data[study] else 1
            if len(totads.split('+')) == 1:
                thisax=ax
        #        if coeff is not None:
        #            print('(%s) %s: %1.2f mV/dec' % (folderlabels[all_folders[study]], study, 1000/coeff[0]))

                if len(this_reddata) > min_datapoints:
                    Tafred.extend([pot_j for pot_j in this_reddata])
                    if a_dat[study]['mode'] == 'COR' and study[:9] != 'this_pcCu':
                        COR_Tafred.extend([pot_j for pot_j in this_reddata])
                    elif study[:9] != 'this_pcCu':
                        CO2R_Tafred.extend([pot_j for pot_j in this_reddata])
                Tafdata.extend([pot_j for pot_j in this_all_data])

                if plot_individual_fits and coeff is not None:
                    thisax.plot(this_reddata[:,0],lin_fun(this_reddata[:,0],*coeff),
                    #plt.plot(this_reddata[:,0],lin_fun(this_reddata[:,0],*coeff),
                            color=studies[study]['color'],
                            markeredgecolor='k',
                            alpha=alpha)

                label='(%s) '%\
                    folderlabels[folder]+'-'.join(study.split('_')[1:])

                if individual_Tafel_slopes_in_label and coeff is not None:
                    label+='%i mV/dec'% (1000/coeff[0])

                #Plot all points containing mass transport limited
                if len(this_reddata):
                    #plt.plot(this_all_data[:, 0], this_all_data[:, 1], markers[all_modes[study]],
                    thisax.plot(this_all_data[:, 0], this_all_data[:, 1], markers[all_modes[study]],
                             markerfacecolor='w', markeredgecolor=studies[study]['color'],
                             alpha=alpha)

                #Plot only mass transport problem free points
                if len(this_reddata):
                    #plt.plot(this_reddata[:, 0], this_reddata[:, 1], markers[all_modes[study]],
                    thisax.plot(this_reddata[:, 0], this_reddata[:, 1], markers[all_modes[study]],
                            markeredgecolor='k' if 'emphasized' in data[study] else None,
                             color=studies[study]['color'],
                             label=label if 'emphasized' in data[study] else None,
                             alpha=alpha,
                             markersize=10 if 'emphasized' in data[study] else 7)
            else:
                tot_all_data[study].extend([pot_j for pot_j in this_all_data])

      if len(totads.split('+')) > 1:
        thisax=ax
        final_all_data, final_red_data = {}, {}
        if len(emphasized_folders):
            nemphstudies=0
            for study in tot_all_data:
                if 'emphasized' in data[study]:
                    nemphstudies+=1
            colors = plt.cm.jet(np.linspace(0, 1, nemphstudies))
            emphcounter=0
        else:
            colors = plt.cm.jet(np.linspace(0, 1, len(tot_all_data.keys())))
        plottedfolders=[]
        for istudy, study in enumerate(tot_all_data.keys()):
            tot_all_data[study] = np.array(tot_all_data[study])
            final_all_data[study] = []

            for pot in np.unique(tot_all_data[study][:, 0]):
              js_at_pot = (np.where(tot_all_data[study][:, 0] == pot)[0])
              final_all_data[study].append([pot, np.log10(np.sum(
                  10**tot_all_data[study][js_at_pot, 1]))])

            final_all_data[study] = np.array(final_all_data[study])

            final_red_data[study], this_all_data, coeff = \
                preselect_datapoints(final_all_data[study],
                                     max_standard_deviation,
                                     potscale,ads,add_pH_to_rate=add_pH_to_rate)

            if all_modes[study] == 'COR':
                        COR_Tafred.extend([pot_j for pot_j in final_red_data[study]])
            else:
                        CO2R_Tafred.extend([pot_j for pot_j in final_red_data[study]])

            if plot_individual_fits and len(final_red_data[study] > 1):
                #plt.plot(final_red_data[study][:,0],lin_fun(final_red_data[study][:,0],*coeff),
                thisax.plot(final_red_data[study][:,0],lin_fun(final_red_data[study][:,0],*coeff),
                        color=studies[study]['color'])

            #plt.plot(final_all_data[study][:, 0],
            thisax.plot(final_all_data[study][:, 0],
                     final_all_data[study][:, 1],markers[all_modes[study]],
                     markerfacecolor='w',
                     markeredgecolor=studies[study]['color'])

            if len(final_red_data[study]) >= min_datapoints:
                #if '-'.join(study.split('-')[:2]) not in plottedfolders:
                if all_folders[study] not in plottedfolders:
                    plottedfolders.append(all_folders[study])
                if verbose:
                    print(string.ascii_uppercase[(len(plottedfolders)-1)],plottedfolders[-1],'%s: %1.2f mV/dec' % (study, 1000/coeff[0]))
                #print(study)
                label='(%s) '%\
                    folderlabels[all_folders[study]]+'-'.join(study.split('_')[1:])
                if individual_Tafel_slopes_in_label:
                    label+='%i mV/dec'% ( 1000/coeff[0])
                all_labels.append(label)
                #plt.plot(final_red_data[study][:, 0],
                thisax.plot(final_red_data[study][:, 0],
                         final_red_data[study][:, 1],markers[all_modes[study]],
                         markeredgecolor='k' if 'emphasized' in data[study] else None,
                         color=studies[study]['color'],
                         label=label if 'emphasized' in data[study] else None,#'%s: %1.1f mV/dec' % (study, 1000/coeff[0]))
                         alpha=alpha,
                        markersize=10 if 'emphasized' in data[study] else 7
                         )
                Tafred.extend([pot_j for pot_j in final_red_data[study]])

            Tafdata.extend([pot_j for pot_j in final_all_data[study]])
      pot_response_fits[totads]=finalize_Tafplot(Tafred, Tafdata, totads, potscale,
              plot_individual_fits, mode,maxpH=maxpH,labels=all_labels,COR=COR_Tafred,
              CO2R=CO2R_Tafred,fit_CO_and_CO2_separately=fit_CO_and_CO2_separately,add_pH_to_rate=add_pH_to_rate,xlim=xlim,ax=thisax,show_rsquared=show_rsquared)
    return studies,folderlabels,all_folders,pot_response_fits


def finalize_Tafplot(Tafred, Tafdata, ads, potscale, plot_individual_fits,mode=None,fit_CO_and_CO2_separately=False,
                     outdir='output', maxpH=100,labels=None,COR=None,CO2R=None,add_pH_to_rate=False,xlim=[-1.69,-0.6],ax=plt,show_rsquared=False):
    thisax=ax
    if fit_CO_and_CO2_separately and COR:
        COR=np.array(COR)
        CORcoeff, dum = curve_fit(lin_fun, COR[:, 0], COR[:, 1])
        thisax.plot(COR[:, 0], lin_fun(COR[:, 0], *CORcoeff), 'b-')
        annotxt='COR: \nTafel slope = %1.1f mV/dec\n'% (1000/CORcoeff[0])
        annotxt+=r'$\alpha$ = %1.2f' % (-CORcoeff[0]/1000*59)
        if show_rsquared:
            annotxt+=', R$^2$=%1.2f'% calculate_rsquared(COR)
#        thisax.annotate('COR: %1.1f mV/dec' % (1000/CORcoeff[0]),
#                 (-1.68,-2.8),fontsize=20).draggable()
        thisax.annotate(annotxt,
                 (-1.62,-3.2),fontsize=20,color='b').draggable()

    if fit_CO_and_CO2_separately and CO2R:
        CO2R=np.array(CO2R)
        CO2Rcoeff, dum = curve_fit(lin_fun, CO2R[:, 0], CO2R[:, 1])
        thisax.plot(CO2R[:, 0], lin_fun(CO2R[:, 0], *CO2Rcoeff), 'k-')
        #thisax.annotate(r' $\alpha$ = %1.2f' % (-CO2Rcoeff[0]*0.059),
                 #(-1.45,0.9),fontsize=20).draggable()
        annotxt='CO2R: \nTafel slope = %1.1f mV/dec\n'% (1000/CO2Rcoeff[0])
        annotxt+=r'$\alpha$ = %1.2f' % (-CO2Rcoeff[0]/1000*59)
        if show_rsquared:
            annotxt+=', R$^2$=%1.2f'% calculate_rsquared(CO2R)
        thisax.annotate(annotxt,(-0.95,1.5),fontsize=20,ha='right',va='top').draggable()
    if '+' not in ads:
            thisax.annotate(ads, (-1.63,-4.4),
                    va='center',ha='left',fontsize=24).draggable()

    if not plot_individual_fits and not fit_CO_and_CO2_separately:
        Tafred, Tafdata = np.array(Tafred), np.array(Tafdata)
        coeff, dum = curve_fit(lin_fun, Tafred[:, 0], Tafred[:, 1])
        print(CRED+'Overall %s Tafel slope: ' % ads, 1000/coeff[0], CEND)
        annotxt='Tafel slope = %1.1f mV/dec\n'% (1000/coeff[0])
        annotxt+=r'$\alpha$ = %1.2f' % (-coeff[0]/1000*59)
        if show_rsquared:
            annotxt+='\nR$^2$=%1.2f'% calculate_rsquared(Tafred)
        thisax.annotate(annotxt,
                 (-0.96,1.6),fontsize=20,ha='right',va='top').draggable()
        thisax.plot(Tafred[:, 0], lin_fun(Tafred[:, 0], *coeff), 'k-')
#        thisax.annotate(r'j$_0$ = 1e%1.2fmA/cm2' % (coeff[0]*(Ueqs[ads]-0.059*7)+coeff[1]),
#                 (-0.65,0.3),fontsize=20,ha='right',va='top').draggable()

    thisax.set_xlabel('U$_{\mathrm{%s}}$ / V' % potscale)


#    if '+' in ads:
#        ylabel='log (j$_{\mathrm{C_{2+}}}$ / mA/cm$^2$)'
#    else:
#        ylabel='log (j$_{\mathrm{CH_4}}$ / mA/cm$^2$)'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outname=outdir+'/Tafel_analysis_%s_maxpH%s_U_%s'% (ads, maxpH, potscale)

    if add_pH_to_rate and ads == 'Methane' and potscale == 'SHE':
            ylabel+=' + pH'
            outname+='+pH'
    if mode:
        outname+='_%s'%mode
    outname+='.pdf'

#    thisax.set_ylabel(ylabel)

#    plt.legend(bbox_to_anchor=(1, 1),prop={'size':6},loc='upper left')
#    plt.legend()
    if not add_pH_to_rate:
        thisax.set_ylim(-5,1.9)
        thisax.set_xlim(xlim)
    #plt.title(ads)
    plt.tight_layout()
    plt.savefig(outname)
    #plt.show()
#    plt.close()
    if not fit_CO_and_CO2_separately:
        return coeff
    elif len(COR):
        return CORcoeff
    elif len(CO2R):
        return CO2Rcoeff

def preselect_datapoints(dat, max_standard_deviation, potscale,ads,add_pH_to_rate=False):
    Tafdata = []
    if isinstance(dat, dict):
        for ipot, pot in enumerate(dat['U_'+potscale]):
            Tafdata.append([pot[0], np.log10(dat['j_partial'][ipot][0])])
        Tafdata = np.array(Tafdata)
    else:
        Tafdata = dat


    if len(Tafdata) < 2:
        return [],[],None

    Tafdata = Tafdata[np.argsort(Tafdata[:, 0])]
    Tafred = detect_mass_trans(Tafdata, max_standard_deviation)
    if add_pH_to_rate and ads == 'Methane' and potscale == 'SHE':
              if potscale == 'SHE':
                Tafred[:, 1] += float(dat['pH'])
                Tafdata[:, 1] += float(dat['pH'])
    if len(Tafred) > 1:
        coeff, dum = curve_fit(lin_fun, Tafred[:, 0], Tafred[:, 1])
    else:
        coeff=None

    return Tafred, Tafdata, coeff


def calculate_rsquared(data):
        coeff, dum = curve_fit(lin_fun, data[:, 0], data[:, 1])
        residuals = data[:, 1] - lin_fun(data[:, 0], *coeff)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((data[:, 1]-np.mean(data[:, 1]))**2)
        r_squared = 1 - (ss_res / ss_tot)
        return r_squared


def detect_mass_trans(Tafdata, max_standard_deviation):
    counter, r_squared = 0, 0
    while r_squared - max_standard_deviation < 0:
        Tafred = Tafdata.copy()
        if counter:
            Tafred = Tafred[counter:]
        if len(Tafred) < 2:
            Tafred=np.array([])
            break
        r_squared = calculate_rsquared(Tafred)
        #coeff, dum = curve_fit(lin_fun, Tafred[:, 0], Tafred[:, 1])
        #residuals = Tafred[:, 1] - lin_fun(Tafred[:, 0], *coeff)
        #ss_res = np.sum(residuals**2)
        #ss_tot = np.sum((Tafred[:, 1]-np.mean(Tafred[:, 1]))**2)
        #r_squared = 1 - (ss_res / ss_tot)
        counter += 1
    return Tafred

def read_data(folders,emphasized_folders):
    data = load_pkl_data(folders)
    if len(emphasized_folders):
        data_emph = load_pkl_data(emphasized_folders)
        for emfol in data_emph.keys():
            data_emph[emfol]['emphasized']=True
        data.update(data_emph)
    else:
        for emfol in data.keys():
            data[emfol]['emphasized']=True

    return data

def lin_fun(x, a, b):
    return a*x+b





if __name__ == "__main__":
    main()
