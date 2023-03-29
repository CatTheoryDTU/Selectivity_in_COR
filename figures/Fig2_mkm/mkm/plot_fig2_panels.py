#!/usr/bin/env python
import sys, os
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import numpy as np
import pickle as pkl

#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
markersize=10

home=os.getcwd()
path_info=home.split('/')
add_total_rate_as_panel=False #Not working yet

nsteps=2

if path_info[-1] == 'complete_run':
    add_total_rate_as_panel=True

if add_total_rate_as_panel: csteps=nsteps+1

fig,ax=plt.subplots(nsteps,2,figsize=(7,7))

def main():
    data,phs,pots=read_data('mkm.pkl')
    R_tot = plot_heatplot(data,phs,pots,nsteps)
    fig.subplots_adjust(wspace=0.05,hspace=0.01)
    fig.tight_layout()
    plt.savefig('Rate_and_selectivities.pdf')
    plot_2D_plot(data,phs,pots,nsteps,potout=[-1.25])

def plot_total_rate(R,ax,X,Y):
    plt.rcParams['figure.figsize'] = (10,5)
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 16
    plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    markersize=10
    fig,ax=plt.subplots()
    b = ax.imshow(R.T,
                interpolation='bicubic',
                cmap=cm.jet,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],norm=LogNorm(),#,
                    vmin=1e-1,
                    vmax=1e5,#)
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())
    fig.colorbar(b,label='Total C$_{2+}$ TOF / s$^{-1}$')
    plt.savefig('Total_rate.pdf')

def plot_2D_plot(data,phs,pots,nsteps,phout=None,potout=None):
    plt.rcParams['figure.figsize'] = (10,5)
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 16
    plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    markersize=10
    fig2d,ax2d = plt.subplots(1,2,figsize=[15,5],sharex=True)

    colors=['b','peru','k']
    lines=['-','--']
    X=np.array(sorted(pots))
    Y=np.array(sorted(phs))

    ax2d[0].plot(np.nan,np.nan,'sb',label='C$_2$H$_4$')
    ax2d[0].plot(np.nan,np.nan,'s',color='peru',label='Oxygenates')
    txtout=''
    if phout is not None:
        xlabel='U$_{\mathrm{SHE}}$ / V'
        #iphs = np.where([i in phout for i in Y])[0]
        rate=np.ones((len(X)))*0.5
        sel=np.ones((len(X)))*0.5
        for iph,ph in enumerate(phout):
          txtout+='# pH %s\n'%ph
          ax2d[0].plot(np.nan,np.nan,'k',label=f'pH{ph}',linestyle=lines[iph])
        #  ax2d[1].plot(np.nan,np.nan,'k',label=f'pH{ph}',linestyle=lines[iph])
          for istep in range(nsteps):
                 for ipot,pot in enumerate(X):
                    rate[ipot]=data[pot][ph][istep]
                    sel[ipot]=rate[ipot]/np.sum(data[pot][ph][:nsteps])

                 ax2d[0].plot(X,rate,color=colors[istep],linestyle=lines[iph],linewidth=2)
                 ax2d[1].plot(X,sel,color=colors[istep],linestyle=lines[iph],linewidth=2)
                 for i,x in enumerate(X):
                     if not istep:
                         txtout+='%1.5f    %1.5f\n'%(x,sel[i])

    elif potout is not None:
        xlabel='pH'
        rate=np.ones((len(Y)))*0.5
        sel=np.ones((len(Y)))*0.5
        for ipot,pot in enumerate(potout):
          txtout+='# pot %s\n'%pot
#          ax2d[0].plot(np.nan,np.nan,'k',label=f'pot {pot}',linestyle=lines[ipot])
        #  ax2d[1].plot(np.nan,np.nan,'k',label=f'pH{ph}',linestyle=lines[iph])
          for istep in range(nsteps):
                 for iph,ph in enumerate(Y):
                    rate[iph]=data[pot][ph][istep]
                    sel[iph]=rate[iph]/np.sum(data[pot][ph][:nsteps])

                 ax2d[0].plot(Y,rate,color=colors[istep],linestyle=lines[ipot],linewidth=2)
                 ax2d[1].plot(Y,sel,color=colors[istep],linestyle=lines[ipot],linewidth=2)
                 for i,x in enumerate(Y):
                     if not istep:
                         txtout+='%1.5f    %1.5f\n'%(x,sel[i])
    else:
        dsds


    ax2d[0].set_ylabel('TOF / s$^{-1}$')
    ax2d[1].set_ylabel('Selectivity')
    ax2d[0].set_xlabel(xlabel)
    ax2d[1].set_xlabel(xlabel)
    ax2d[0].set_yscale('log')
#    ax2d[0].set_ylim([1e-4,10000])
    ax2d[1].set_ylim([0,1])
#    ax2d[0].set_xlim([-1.6,-1.0])
#    ax2d[0].legend()
#    ax2d[1].legend()
    fig2d.tight_layout()
    plt.savefig('2D_selectivities.pdf')
#    plt.show()

    out=open('Selectivities.out','w')
    out.write(txtout)
    out.close()



def plot_heatplot(data,phs,pots,nsteps):
    X=np.array(sorted(pots))
    Y=np.array(sorted(phs))

    R_tot=np.ones((len(X),len(Y)))
    for col in range(2):
     for istep in range(nsteps):
        R,S = get_rate_and_selectivity(col,istep,data,nsteps,X,Y)
        plot_it(R,S,ax,col,istep,X,Y)
        if istep == nsteps-1:
            for thisax in ax[istep]:
                thisax.set_xlabel('U$_{\mathrm{SHE}}$ / V')
        else:
            for thisax in ax[istep]:
                thisax.set_xticks([])

        ax[istep][0].set_ylabel('pH')
        ax[istep][1].set_yticks([])
        R_tot+=R
    return R_tot




def get_rate_and_selectivity(col,istep,data,steps,X,Y):

    Selectivity=np.ones((len(X),len(Y)))*0.5
    rate=np.ones((len(X),len(Y)))*0.5

    for ix,x in enumerate(X):
       for iy,y in enumerate(Y):
        try:
            if col == 1:
                Selectivity[ix][iy]=data[x][y][istep]/np.sum(data[x][y][:nsteps])
            else:
                rate[ix][iy]=data[x][y][istep]#/np.sum(data[x][y][:nsteps])
        except:
            Selectivity[ix][iy]=np.nan#data[x][y][istep]#/np.sum(data[x][y][:3])
            rate[ix][iy]=1e-20#np.nan#data[x][y][istep]#/np.sum(data[x][y][:nsteps])
            #print(x,y)
            pass
    return rate, Selectivity

def plot_it(R,S,ax,col,istep,X,Y):
        thisax=ax
        if len(ax) >1: thisax=ax[istep][col]
        if col == 0:
         b = thisax.imshow(R.T,
                interpolation='bicubic',
                cmap=cm.jet,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],norm=LogNorm(),#,
                    vmin=1e-3,
                    vmax=1e5,#)
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())

        else:
         a = thisax.imshow(S.T,
                interpolation='bicubic',
                cmap=cm.RdYlGn,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],#norm=LogNorm(),#,
                    vmin=0,
                    vmax=1,
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())


def read_data(infile='mkm.pkl'):
    data_in = pkl.load(open(infile,'rb'),encoding='latin1')
    data={}
    pots,phs=[],[]
    for dat in data_in['rate_map']:
        pot,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot not in data:
            data[pot] = {}
        data[pot][ph] = dat[1]
        if pot not in pots:
            pots.append(pot)
        if ph not in phs:
            phs.append(ph)
    return data,phs,pots
if __name__ == "__main__":
    main()
