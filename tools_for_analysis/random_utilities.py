import os
from ase.io import read,write
import numpy as np

def write_pckl_data(path,data):
        import pickle
        ff=open(path,'wb')
        pickle.dump(data, ff)
        ff.close()

def create_new_dir(home,name):
    if name not  in os.listdir(home):
        os.mkdir(name)
    else:
        x = input('Structure directory already exists. Continue?')
        if x in ['n','N','No','q']:
            print('Neglecting it')
            pass

def create_dir_and_save_old_output(outdir,move_wavecar=False,move_contcar=False):
    home = os.getcwd()
    latest_run=None

    if outdir[-1] == '/':
        outdir = outdir[:-1]

    if outdir in os.listdir(home):
        #Check if an old run is already there
        for i in range(50):
            if outdir+'_%i'%i in os.listdir(home):
                pass
            else:
                latest_run=i
                break
    if latest_run is not None:
        os.system('mv   %s   %s'%(outdir,outdir+'_%i'%latest_run))

    os.mkdir(outdir)

    if latest_run is not None and move_contcar:
        try:
            lastat=read(outdir+'_%i/vasprun.xml'%latest_run)
        except IndexError:
            pass
        else:
            initat=read('init.traj')
            posdif = lastat.get_positions() - initat.get_positions()
            if np.sum(np.linalg.norm(posdif,axis=1)) > 0.2:
                os.system('mv init.traj init_bakk.traj')
                write('init.traj', lastat)

    if move_wavecar:
        if latest_run is not None:
            os.system('mv %s/WAVECAR %s'%(outdir+'_%i'%latest_run,outdir))
            os.system('mv %s/CHGCAR %s'%(outdir+'_%i'%latest_run,outdir))
            os.system('mv %s/CHG %s'%(outdir+'_%i'%latest_run,outdir))
            os.system('mv %s/DOSCAR %s'%(outdir+'_%i'%latest_run,outdir))

def get_reference_vibrational_contribution(adsorbate,CO_as_ref=True,code='GPAW',references={'C':'CO'}):
    E={}
    if code == 'GPAW':
        E['O'] = get_gas_reference_vibrational_contributions('H2O')-get_gas_reference_vibrational_contributions('H2')
        E['H'] = 0.5*get_gas_reference_vibrational_contributions('H2')
        E['CO'] = get_gas_reference_vibrational_contributions('CO')
        #E['C'] = get_gas_reference_vibrational_contributions('CO')-E['O']
        Cref=''
        for ilet, letter in enumerate(references['C']):
            if letter.isnumeric():
                if ilet == 0:
                    print ('C reference starts with a number, can not be interpreted')
                    return False
                for nat in range(int(letter)-1):
                    Cref+=references['C'][ilet-1]
            else:
                Cref+=letter
        E['C']=get_gas_reference_vibrational_contributions(references['C'],code='GPAW')
        for letter in Cref:
            if letter == 'C': continue
            E['C']-=E[letter]
    else:
        print('Only GPAW is implemented for vibrational contribution of gas species at the moment')
        sys.exit()

    #E_OCCO = 2*(E_CO)
    E_ref=0
    #print('ads',adsorbate)
    for elem in adsorbate:
        E_ref += E[elem]

    return E_ref

def get_gas_reference_vibrational_contributions(gas,code='GPAW'):
    if code == 'GPAW':
        gases={'CH4': 0.6609615722468455, #1atm
                'H2':-0.05276775930248084, #1atm
                'H2O': 0.0071360210575864835, #0.035atm
                'CO': -0.38954084777493236, #1atm
                'CO2':-0.25199999999999534 #1atm
                }
    elif code == 'VASP':
        gases={'CH4': 0.6609615722468455, #1atm
                'H2':-0.035, #1atm
                'H2O': 0.009, #0.035atm
                'CO': -0.371 #1atm
                }
    else:
        raise InputError('Code not understood')
    return gases[gas]

def get_reference_energies(adsorbate,CO_as_ref=True,code='VASP',use_potential_energy=True,
        references={'C':'CO'}):
    E = {}
    #OCCO (g)
    if code == 'QE':
        E['C'] = get_gas_reference('CH4')-2*get_gas_reference('H2')
        E['O'] = get_gas_reference('H2O')-get_gas_reference('H2')
        E['H'] = 0.5*get_gas_reference('H2')
        E['CO'] = get_gas_reference('CO')
    elif code == 'VASP':
        E['H'] = 0.5*get_gas_reference_VASP('H2')
        E['O'] = get_gas_reference_VASP('H2O')-get_gas_reference_VASP('H2')
        Cref=''
        for ilet, letter in enumerate(references['C']):
            if letter.isnumeric():
                if ilet == 0:
                    print ('C reference starts with a number, can not be interpreted')
                    return False
                for nat in range(int(letter)-1):
                    Cref+=references['C'][ilet-1]
            else:
                Cref+=letter
        E['C']=get_gas_reference_VASP(references['C'])
        for letter in Cref:
            if letter == 'C': continue
            E['C']-=E[letter]

        #E['C'] = get_gas_reference_VASP('CH4')-2*get_gas_reference_VASP('H2')
        #E['O'] = get_gas_reference_VASP('H2O')-get_gas_reference_VASP('H2')
        #E['H'] = 0.5*get_gas_reference_VASP('H2')
        #E['CO'] = get_gas_reference_VASP('CO')
        #if CO_as_ref:
        #    E['C'] = get_gas_reference_VASP('CO')-E['O']
    elif code == 'GPAW':
        E['O'] = get_gas_reference_GPAW('H2O',use_potential_energy)-get_gas_reference_GPAW('H2',use_potential_energy)
        E['H'] = 0.5*get_gas_reference_GPAW('H2',use_potential_energy)
        Cref=''
        for ilet, letter in enumerate(references['C']):
            if letter.isnumeric():
                if ilet == 0:
                    print ('C reference starts with a number, can not be interpreted')
                    return False
                for nat in range(int(letter)-1):
                    Cref+=references['C'][ilet-1]
            else:
                Cref+=letter
        E['C']=get_gas_reference_GPAW(references['C'])
        for letter in Cref:
            if letter == 'C': continue
            E['C']-=E[letter]
        E['CO'] = get_gas_reference_GPAW('CO',use_potential_energy)

    #E_OCCO = 2*(E_CO)
    E_ref=0
    #print('ads',adsorbate)
    for elem in adsorbate:
        E_ref += E[elem]
    #print(adsorbate,E_ref)

    return E_ref

def add_magnetic_moment(atoms):
    melements= ['Ni','Fe','Mo','Mn','Ti','V','Cr','Co']
    for atom in atoms:
        if atom.symbol == 'Ti':
                atom.magmom = 3.
        elif atom.symbol in melements:
                atom.magmom = 5.
        else:
                atom.magmom = 0.

    if np.any(np.around(atoms.arrays['initial_magmoms'],2) > 0.):
        spinpol = True
    else:
        spinpol = False
    return atoms,spinpol

def metal_reference_energy(element):
       ene = {'Ag':  0.55065969,
              'Au':  -0.26473304,
              'Cu':  -0.64186417,
              'Ir':  -5.7997582,
              'Ni':  -2.53057222,
              'Pd':  -1.94717605,
              'Pt':  -3.1320257,
              'Rh':  -4.24356899,
              'Zn':  2.0471165,
              'Sn': -1.20811615,
              'Tl': 0.9362309,
              'Pb': -.83958187,
              'Bi': -1.1706735,
              'Sb': -2.08191495,
              'Hg': 2.4382504,
              'Ga': -1.5061499,
              'Ge': -3.05828215,
              'In': -.71689777,
              'Cd': 2.51238375,
              'Hf': -7.947293,
              'Ti': -5.888111,
              'Sc': -3.6235767,
              'Zr': -5.824230885,
              'Y': -4.02269914,
              'W': -5.2008225,
              'Fe': -5.64679155,
              'Cr': -7.175663,
              'Mn': -6.407885810689655,
              'Co': -4.294195545,
              'V': -6.62063308,
              'Ta': -9.39766076,
              'Mo': -8.14826917,
              'Nb': -7.25114749,


            }
       return ene[element]

def metal_surface_energies(element,facet):
    surf = {'Ag': {'100': 0.051, '111': 0.048,'211': 0.054},
            'Au': {'100': 0.054, '111': 0.046, '211': 0.051},
            'Cu': {'100': 0.092, '111': 0.082, '211': 0.102},
            'Ir': {'100': 0.176, '111': 0.143, '211': 0.169},
            'Ni': {'100': 0.138, '111': 0.120, '211': 0.140},
            'Pd': {'100': 0.095, '111': 0.084, '211': 0.101},
            'Rh': {'100': 0.145, '111': 0.124, '211': 0.144},
            'Pt': {'100': 0.115, '111': 0.092, '211': 0.110},
            'Zn': {'111': 0.022, '100':0.033},
            'Cd': {'111': 0.012, '100':0.017},
            'Sn': {'100': 0.055, '111': 0.041},
            'Hg': {'100': 0.004, '111': 0.004},

            }
    # Actually it's: 'Zn': {'0001': 0.022, '10-10':0.033},
    #Careful: Sn (diamond structure has much more stable terminations!)
    return surf[element][facet]

#Prerelaxation for vaspsol
def prerelax(atoms,home,mirror=False,fmax=0.1,k=3,unstable=False,constraint='internal',dmax=3):
    Cind,Oind=[],[]
    for atmind,atom in enumerate(atoms):
                if atom.symbol == 'C':
                    Cind.append(atmind)
                elif atom.symbol == 'O':
                    Oind.append(atmind)
    nOCCO,nCO2=0,0
    if unstable:
        if len(Cind) == len(Oind): #For OCCO
            if len(Cind) > 1:
                for i,Cindex in enumerate(Cind):
                 for j,Cindex2 in enumerate(Cind):
                  if Cindex != Cindex2:
                    C1pos=atoms.positions[Cindex]
                    C2pos=atoms.positions[Cindex2]

                    zdist=np.around(abs(C1pos[2]-C2pos[2]),3)
                    xycell = np.diag(atoms.cell.copy())
                    if zdist < dmax:#zdist > 0.0:
                        for nx in [-1,0,1]:
                           for ny in [-1,0,1]:
                                xydist = abs(np.linalg.norm(C1pos-(C2pos+np.array([nx,ny,0])*xycell)))
                                        #print(zdist,nx, ny, xydist)
                                if xydist < dmax:
                                    nOCCO += 1
                                    dist = np.linalg.norm(C1pos-(C2pos+np.array([nx,ny,0])*xycell))
                nOCCO/=2
            else:
                dist=100.

            if len(Cind) > 2 and not mirror:
                if np.around(dist - np.linalg.norm(atoms.positions[Cind[2]]-atoms.positions[Cind[3]]),3) == 0. and\
                    atoms.positions[Cind[0]][2] - atoms.positions[Cind[2]][2] > 5:
                        raise Exception('It seems like the structure is symmetric, but mirror is set to False. Please check!')

        elif 2*len(Cind) == len(Oind): #CO2
            if len(Cind):
                angle_const_pairs=[]
                nearest_mets=[]
                CMdists=[]
                Oindex_sorted_by_C=[]
                Odists=[]
                for iC,Cindex in enumerate(Cind):
                    Cpos=atoms.get_positions()[Cindex]
                    Oneighs=[]
                    for Oindex in Oind:
                        Opos = atoms.get_positions()[Oindex]-Cpos
                    #Opos=np.array([atoms.get_positions[Oindex]-Cpos for Oindex in Oind])

                        zdist=abs(Opos[2])
                        if zdist > 3: continue
                        else:
                            Oneighs.append(Oindex)

                    if len(Oneighs) != 2:
                        raise Exception('There do not seem 2Os being connected to a carbon')
                    nCO2+=1
                    Oindex_sorted_by_C.append([Oneighs[0],Oneighs[1]])
                    Odists.append(np.linalg.norm(atoms.get_positions()[Oneighs[1]]-atoms.get_positions()[Oneighs[0]]))

                    angle_const_pairs.append([Oneighs[0],Cindex,Oneighs[1]])

                    # Find nearest metal atom
                    dists=atoms.get_positions()-Cpos
                    sortdist=np.argsort(np.linalg.norm(dists,axis=1))
                    for dist in sortdist:
                        if atoms.get_chemical_symbols()[dist]  not in ['C','O']:
                            nearest_mets.append(dist)
                            CMdists.append(np.linalg.norm(dists[dist]))
                            break

            else:
                dist=100.
                nCO2=0

    else:
        dist=100
        nOCCO=0

    if any([nCO2,nOCCO]) and unstable:
     if nCO2 and nOCCO:
         raise Exception('It appears both OCCO and CO2 are in teh same cell! Aborted')
     if nCO2:
         #Still experimental at this point
         if constraint == 'internal':
            from ase.constraints import FixInternals
            from math import pi
            CMbond = [CMdists[0], [Cind[0], nearest_mets[0]]]
            OObond = [Odists[0],Oindex_sorted_by_C[0]]
            angle = [atoms.get_angle(*angle_const_pairs[0]), angle_const_pairs[0]]
            if mirror:
                angle2 = [atoms.get_angle(*angle_const_pairs[1]), angle_const_pairs[1]]
                CMbond2 = [CMdists[1], [Cind[1], nearest_mets[1]]]
                OObond2 = [Odists[1],Oindex_sorted_by_C[1]]

                #c = FixInternals(bonds=[CMbond,CMbond2,OObond,OObond2],
                c = FixInternals(bonds=[OObond,OObond2],
#                     angles=[angle,angle2])
                        )
            else:
                #c = FixInternals(bonds=[CMbond,OObond],# angles=[angle1],
                c = FixInternals(bonds=[OObond],# angles=[angle1],
                     #angles=[angles])
                    )
     if nOCCO:
        if constraint == 'internal':
            from ase.constraints import FixInternals
            from math import pi
            bondC1 = [np.linalg.norm(atoms.positions[Cind[0]]-atoms.positions[Cind[1]]), [Cind[0], Cind[1]]]
            bondO1 = [np.linalg.norm(atoms.positions[Oind[0]]-atoms.positions[Oind[1]]), [Oind[0], Oind[1]]]
            dihedral_indices1 = [Oind[0],Cind[0], Cind[1], Oind[1]]
            dihedral1 = [atoms.get_dihedral(*dihedral_indices1) * pi / 180,
                            dihedral_indices1]
            if mirror:
                bondC2 = [np.linalg.norm(atoms.positions[Cind[0]]-atoms.positions[Cind[1]]), [Cind[2], Cind[3]]]
                bondO2 = [np.linalg.norm(atoms.positions[Oind[0]]-atoms.positions[Oind[1]]), [Oind[2], Oind[3]]]
                dihedral_indices2 = [Oind[2],Cind[2], Cind[3], Oind[3]]
                dihedral2 = [atoms.get_dihedral(*dihedral_indices1) * pi / 180,
                            dihedral_indices2]
                c = FixInternals(bonds=[bondC1,bondC2,bondO1,bondO2],# angles=[angle1],
                     dihedrals=[dihedral1,dihedral2])

            else:
                c = FixInternals(bonds=[bondC1,bondO1],# angles=[angle1],
                     dihedrals=[dihedral1])

        elif constraint == 'bondlength':
            from ase.constraints import FixBondLengths
            if mirror:
                c=FixBondLengths([[Cind[0],Cind[1]],[Oind[0],Oind[1]],
                                  [Cind[0],Oind[1]],[Cind[1],Oind[0]],
                                  [Cind[2],Cind[3]],[Oind[2],Oind[3]],
                                  [Cind[2],Oind[3]],[Cind[3],Oind[2]]])

            else:
                c=FixBondLengths([[Cind[0],Cind[1]],[Oind[0],Oind[1]],
                                  [Cind[0],Oind[1]],[Cind[1],Oind[0]]])

        elif constraint == 'spring':
            #Don't use it
            from ase.constraints import Hookean
            c = Hookean(a1=Cind[0],a2=Cind[1],rt=dist, k=k)
 #           cO = Hookean(a1=Cind[0],a2=Cind[1],rt=dist, k=k)
            nconst=1
        else:
            raise TypeError('Prerelax constraint type not understood')
     atoms.constraints.append(c)

    outdir = 'prerelax'
    create_dir_and_save_old_output(outdir)
    os.system('cp run.py %s'%outdir)
    os.chdir(outdir)
    atoms.calc.set(lsol=False)
    if any([nOCCO,nCO2]) and unstable:
        print(atoms.constraints)
        from ase.optimize import BFGS
        nsw_bakk = atoms.calc.int_params['nsw']
        atoms.calc.set(nsw=0)
        opt=BFGS(atoms,logfile='relax.log',trajectory='relax.traj')
        opt.run(fmax=fmax)
        del atoms.constraints[-1]
        #for dummy in range(nconst):
        #    del atoms.constraints[-1]
        #    if mirror:
        #        del atoms.constraints[-1]
        atoms.calc.set(nsw=nsw_bakk)
    else:
        ediffg_bakk = atoms.calc.exp_params['ediffg']
        atoms.calc.set(ediffg=-0.1)
        energy = atoms.get_potential_energy()
        atoms.calc.set(ediffg=ediffg_bakk)
    os.system('rm CHG')
    atoms.calc.set(lsol=True,lambda_d_k=3.0,tau=0)
    os.chdir(home)
    #write('Prerelax.traj',atoms)
    return atoms


def check_if_OCCO_split(atoms,dmax=2.5):
                C_pos = {}
                for i,pos in enumerate(atoms.positions):
                    if atoms.get_chemical_symbols()[i] == 'C':
                        C_pos[i]=pos

                nOCCO=0
                for i,Cind in enumerate(C_pos.keys()):
                    for j,C2ind in enumerate(C_pos.keys()):
                        C,C2=C_pos[Cind],C_pos[C2ind]
                        zdist=np.around(abs(C[2]-C2[2]),3)
                        if zdist < dmax and Cind != C2ind:#zdist > 0.0:
                            for nx in [-1,0,1]:
                             for ny in [-1,0,1]:
                                xycell = np.diag(atoms.cell.copy())

                                xydist = np.linalg.norm(C-(C2+np.array([nx,ny,0])*xycell))
                                #print(zdist,nx, ny, xydist)
                                if xydist < dmax:
                                    nOCCO += 1

                return nOCCO

def check_if_CO2desorbed(atoms,dmax=3.0,mirror=True):
                C_pos = {}
                for i,pos in enumerate(atoms.positions):
                    if atoms.get_chemical_symbols()[i] == 'C':
                        C_pos[i]=pos

                nCO2=0
                for i,Cind in enumerate(C_pos.keys()):
                    allpos=[]
                    for x in [-1,0,1]:
                        for y in [-1,0,1]:
                            dist=np.linalg.norm(atoms.get_positions()+np.diag(atoms.get_cell())*np.array([x,y,0])-C_pos[Cind],axis=1)
                            sortdist=np.argsort(dist)
                            for ineigh in sortdist:
                                if atoms.get_chemical_symbols()[ineigh] in ['C','O']:
                                    continue
                                zdist = np.around(abs(dist[ineigh]),3)
                                if zdist < dmax:
                                    nCO2+=1
                                    print(zdist,ineigh)
#                            break
                        #else:
                        #Put the following in if you get problems with pbc
#                            for nx in [-1,0,1]:
#                             for ny in [-1,0,1]:
#                                xycell = np.diag(atoms.cell.copy())
#                                xydist = np.linalg.norm(atoms.get_positions()[ineigh]-(C_pos[i])*xycell))
#
#                                xydist = np.linalg.norm(C-(C2+np.array([nx,ny,0])*xycell))
#                                #print(zdist,nx, ny, xydist)
#                                if xydist < dmax:
#                                    nOCCO += 1
                if mirror:
                    return nCO2//2
                else:
                    return nCO2

def detect_surface_coordination(image,zval):
            surfpos=[]
            for posind,pos in enumerate(image.positions):
                if abs(pos[2] - zval[-1]) < 1e-3:
                        surfpos.append(pos[:2])
            surfpos-=np.linalg.norm(image.cell,axis=1)[:2]*0.5

            surfpos=np.array(surfpos)-surfpos[np.argsort(np.linalg.norm(surfpos,axis=1))][0]
            surfpos = surfpos[np.argsort(np.linalg.norm(surfpos,axis=1))]
            v1,v2 = surfpos[1],surfpos[2]
            angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            if np.around(angle/(2*np.pi)*360,0) == 180:
                v1,v2 = surfpos[1],surfpos[3]
                angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            return angle/(2*np.pi)*360,[np.linalg.norm(v1),np.linalg.norm(v2)]
            #continue

def get_crystal_from_info(dire='..'):
    info = open(dire+'/info','r')
    infoline=info.readlines()[2].split()
    for i,inf in enumerate(infoline):
        if inf == "'crystal_system':":
            lattice = infoline[i+1].lstrip("'").rstrip("',")
        if inf == "'symbol':":
            pointgroup = infoline[i+1].lstrip("'").rstrip("',")

    return lattice,pointgroup

def get_gas_reference_QE(gas):
    gases = { 'CH4': -231.5729033119903,
              'H2': -32.941895798946696,
              'H2O': -496.2715343769096,
              'CO': -626.4870639805104
              }
    return gases[gas]

def get_gas_reference_GPAW(gas,use_potential_energy=True):
    if use_potential_energy:
        gases = { 'CH4': -34.082626231,
                  'H2': -8.009969,
                  'H2O': -27.231564,
                 'CO': -34.797848,
                 'CO2': -54.767 #+ 0.33 Empirical correction only for OCO backbone
                  }
    else:
        gases={'CH4': -33.422,
                'H2':-8.063,
                'H2O':-27.224,
                'CO': -35.187,
                'CO2': -55.019  #+0.33#Empirical correction for OCO backbone
                }
    return gases[gas]

def get_gas_reference_VASP(gas):
    #All are potential energies
    gases = { 'CH4': -23.279576,
              'H2': -7.1676765,
              'H2O': -12.806973,
              'CO': -12.067233,
              'CO2': -18.401744 + 0.45
              }
    return gases[gas]

def lin_fun(x,a,b):
    return a*x+b

def quad_fun(x,a,b,c):
    return a*x*x+b*x+c

def get_reduced_alloy_name(comp):
            #print(comp)
            if len(comp) < 3:
                return comp
            #overlayers
            if '_on_' in comp: return comp
            if comp in ['Au4Cu8Pd4']:
                return 'AuCu2Pd'

            compl=list(comp)
            for il,letter in enumerate(compl):
                try:
                    dummy=int(letter)
                except ValueError:
                    compl[il] = ' '
            letterinds=np.where(np.array(compl) == ' ')[0]
            stoich = np.array([int(i) for i in ''.join(compl).split()])

            #print(compl)
            while not np.any(stoich%2):
                        stoich = stoich/2#0.5

            if np.all(stoich > 1):
                    if stoich[1] > stoich[0] and not stoich[1]%stoich[0]:
                        stoich//=stoich[0]
                    elif stoich[0] >= stoich[1] and not stoich[0]%stoich[1]:
                        stoich//=stoich[1]

            stoich=list(stoich)
            for nati,nat in enumerate(stoich):
                  if int(nat) == 1:
                      stoich[nati] = ''
                  else:
                      stoich[nati] = int(stoich[nati])
            outname=''
            istoich=0
            for i in range(10):
                if i in letterinds:
                    outname+=comp[i]
                else:
                    outname+=str(stoich[istoich])
                    istoich+=1
                if istoich == len(stoich):
                    break

            #print(comp,outname)
            #das
            return outname#'%s%s%s%s'%(comp[:2],str(stoich[0]),comp[3:5],str(stoich[1]))

            #old routine delete in case the new one works
            if comp not in ['Au4Cu8Pd4'] and len(comp) > 2:
                if len(comp) == 6:
                    stoich = np.array([float(comp[2]),float(comp[5])])
                    stoich = np.array([int(comp[2]),int(comp[5])])
                else:
                    try:
                        stoich = np.array([float(comp[2]),float(comp[-2:])])
                        stoich = np.array([int(comp[2]),int(comp[-2:])])
                    except:
                        stoich = np.array([float(comp[2:4]),float(comp[-1])])
                        stoich = np.array([int(comp[2:4]),int(comp[-1])])

                while not np.any(stoich%2):
                        stoich = stoich/2#0.5

                if np.all(stoich > 1):
                    if stoich[1] > stoich[0] and not stoich[1]%stoich[0]:
                        stoich//=stoich[0]
                    elif stoich[0] >= stoich[1] and not stoich[0]%stoich[1]:
                        stoich//=stoich[1]

                stoich=list(stoich)
                for nati,nat in enumerate(stoich):
                  if int(nat) == 1:
                      stoich[nati] = ''
                  else:
                      stoich[nati] = int(stoich[nati])
            else:
                if len(comp) < 3:
                    return comp
                else:
                    return 'AuCu2Pd'
            return '%s%s%s%s'%(comp[:2],str(stoich[0]),comp[3:5],str(stoich[1]))

def update_dicts(alldicts,data):
    for dicti in alldicts:
        dicti[data]={}

def get_cellsize(atoms):
    return np.product(np.diag(atoms.cell[:2,:2].copy()))*1e-16

def get_bands_DOS(pdos):
      if len(pdos) == 10: # spin unpolarised
          s, p_y, p_z, p_x, d_xy, d_yz, d_z2, d_xz, d_x2_y2 =  pdos[1:,1:]
          data = {'s':s, 'p_y':p_y, 'p_z':p_z, 'p_x':p_x, 'd_xy':d_xy, 'd_yz':d_yz,\
                  'd_z2':d_z2, 'd_xz':d_xz, 'd_x2_y2':d_x2_y2}
          return data

      elif len(pdos) == 19: # spin polarized
          s_up, s_down, p_y_up, p_y_down, p_z_up, p_z_down, p_x_up, p_x_down, \
                  d_xy_up, d_xy_down, d_yz_up, d_yz_down, d_z2_up, d_z2_down, \
                  d_xz_up, d_xz_down, d_x2_y2_up, d_x2_y2_down =  pdos[1:,1:]
          data = {'s_up':s_up, 's_down':s_down, 'p_y_up':p_y_up, 'p_y_down':p_y_down,\
                  'p_z_up':p_z_up, 'p_z_down':p_z_down, 'p_x_up':p_x_up, 'p_x_down':p_x_down,\
                  'd_xy_up':d_xy_up, 'd_xy_down':d_xy_down, 'd_yz_up':d_yz_up, 'd_yz_down':d_yz_down, \
                  'd_z2_up':d_z2_up, 'd_z2_down':d_z2_down, 'd_x2_y2_up':d_x2_y2_up, 'd_x2_y2_down':d_x2_y2_down,\
                  'd_xz_up':d_xz_up, 'd_xz_down':d_xz_down}
          return data
      print("Missing data from DOSCAR - check calculation")

def d_band_info(energies, pdos):
    from scipy.integrate import quad, simps
    def integrand_numerator(e, rho):
        return e * rho
    def integrand_denominator(e,rho):
        return rho

    # Getting mean d-band
    epsilon_d_num = simps( integrand_numerator(energies,pdos), energies)
    epsilon_d_den = simps(integrand_denominator(energies,pdos), energies)
    epsilon_d = epsilon_d_num / epsilon_d_den

    # First moment of d-band
    def integrand_moment_numerator(e, epsilon_d, moment, rho):
        return ( (e - epsilon_d ) ** moment ) * rho
    w_d_numerator = simps( integrand_moment_numerator(energies, epsilon_d, 2, \
            pdos), energies)
    w_d_denominator = epsilon_d_num
    w_d = w_d_numerator / w_d_denominator

    return [epsilon_d, w_d]

def check_stoichiometry_alloys(atoms,infofile):
    info = open(infofile,'r')
    infolines=info.readlines()
    comp = infolines[0]

    stoichlist=[]
    for i,s in enumerate(comp):
        if s.isdigit():
            stoichlist.append(int(s))
    #    else:
    #        stoichlist.append(0)
    if len(stoichlist) > 2:
        if comp[-3].isdigit():
            stoichlist=[stoichlist[0],stoichlist[1]*10+stoichlist[2]]
        else:
            stoichlist=[stoichlist[0]*10+stoichlist[1],stoichlist[2]]

#    if len(stoichlist) > 2:
#        print('Stoichiometry error!')
#        sys.exit()

    element_ratio = np.unique(atoms.get_atomic_numbers(),return_counts=True)[1]
    stoichiometry = stoichlist[0]/np.sum(stoichlist)

    if (np.around(element_ratio[0]/np.sum(element_ratio),3) == np.around(stoichiometry,3) or\
        np.around(1 - element_ratio[0]/np.sum(element_ratio),3) == np.around(stoichiometry,3)):\
        return True,stoichlist
    else:
        return False,stoichlist
