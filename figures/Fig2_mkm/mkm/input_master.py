import sys

# General settings - usually unchanged
scaler = 'ThermodynamicScaler'
descriptor_names= ['voltage', 'pH']
#descriptor_ranges= [[-1.5,0],[3,8]]
descriptor_ranges= [[-1.7,-0.9],[7,14]]
temperature = 300.
#resolution=[15,15]
resolution=[17,15]

gas_thermo_mode = 'frozen_gas'
adsorbate_thermo_mode = 'frozen_adsorbate'
#electrochemical_thermo_mode = "simple_electrochemical"
electrochemical_thermo_mode = 'surface_charge_density'

# solver settings
decimal_precision = 100
tolerance = 1e-25
max_rootfinding_iterations = 200
max_bisections = 3

interaction_strength=float(sys.argv[1])

beta = 0.5
# Cu - CO2 reduction input file settings
#input_file = '../energy_input.txt'
input_file = '../../../parsed_data_for_analysis/mkm_energy_input_freEn.txt'
surface_names = ['Cu']
potential_reference_scale = 'SHE'

species_definitions = {}
# pressures
#species_definitions['H_g'] = {'pressure':1.0}
species_definitions['ele_g'] = {'pressure':1.0, 'composition':{}}
species_definitions['CO_g'] = {'pressure':1.0}
species_definitions['H2_g'] = {'pressure':0.0}
species_definitions['H2O_g'] = {'pressure':0.035}
species_definitions['CH4_g'] = {'pressure':0.0}
species_definitions['O2_g'] = {'pressure':0.0}
species_definitions['CH3CH2OH_g'] = {'pressure':0.0}
species_definitions['CH3COOH_g'] = {'pressure':0.0}
species_definitions['CH3COO_g'] = {'pressure':0.0}
species_definitions['C2H4_g'] = {'pressure':0.0}
species_definitions['OH_g'] = {'pressure':0.0}

estimate_frequencies = 0

# interaction model

adsorbate_interaction_model = 'first_order' #use "single site" interaction model
#adsorbate_interaction_model = None #use "single site" interaction model

interaction_response_function = 'smooth_piecewise_linear' #use "smooth piecewise linear" interactions

interaction_fitting_mode = None
cross_interaction_mode = 'geometric_mean' #use geometric mean for cross parameters

interaction_scaling_constraint_dict = {
                                'CO_211':[0,0,None],
                                'H_211':[0,0,None],
                                'OH_211':[0,0,None],
                                'COH_211':[0,0,None],
                                'CHO_211':[0,0,None],
                                'CH_211':[0,0,None],
                                'CH2_211':[0,0,None],
                                'CH3_211':[0,0,None],
                          }

if 1 == len(surface_names):
        numBeforePt = 0
        numAfterPt = 0
else:
        numBeforePt = len(surface_names)-surface_names.index('Pt')
        numAfterPt = len(surface_names)-numBeforePt-1
        transition_state_cross_interaction_mode = 'transition_state_scaling' #use TS scaling for TS interaction

# site interaction scaling:
eHCO=0.7274
eCO=2.4670*interaction_strength
eH=0.0

species_definitions['100'] = {'site_names': ['100'], 'total':1.0}
species_definitions['100']['interaction_response_parameters'] = {'cutoff':0.25,'smoothing':0.05}
species_definitions['dl'] = {'site_names': ['dl'], 'total':1.0}
species_definitions['CO_100'] = {
                'self_interaction_parameter':[None]*numBeforePt+[eCO]+[None]*numAfterPt,
}

species_definitions['OCCO_100'] = {
                   'self_interaction_parameter':[None]*numBeforePt+[1*eCO]+[None]*numAfterPt,
                   'n_sites':2
                   }


for sp in ['CO_100','CH2OH_100','O_100','OH_100',]:
    species_definitions[sp] = species_definitions['CO_100'].copy()

for sp in ['CO-CO_100','COCO_100','OCCO_100','OCCO-H2O-ele_100','OCCOH_100','OCCOH-ele_100','OCCOH-H2O-ele_100','HOCCOH_100','HOCCOH-ele_100','CCO_100','CCO-H2O-ele_100','HOCCHO_100','HOCCHOH_100','HOCCH2O_100','OCCH2OH_100','CCOH_100','CCOH-ele_100','HCCO_100','HCCO-H2O-ele_100','CCHO_100','CCH2O_100','CC_100','CC-H2O-ele_100','H2CCO_100','H2CCO-H2O-ele_100','HCCOH_100','HCCOH-ele_100','HCCHO_100','CCH_100','OHCCH2_100','H2CCHO_100','HCCHOH_100','H2CCOH_100','HCCH_100','OCHCH3_100','H2CCH2O_100','H2CCH_100','H3CCH2O_100','OCC-H2O-ele_100','OHCC-H2O-ele_100','CHCO-H2O-ele_100','CH2CO-H2O-ele_100','H3CCO_100','CCH2O-H2O-ele_100','H2CCH2O-H2O-ele_100','OH2CC-H2O-ele_100','COC-H2O-ele_100','COHC-H2O-ele_100', 'H2CCO-H2O-ele_100', 'OH2CC-H2O-ele_100', 'COH2C-H2O-ele_100','CCH2_100','HCC-OH-ele_100','HCCOH-H2O-ele_100','HOCCH3_100','CCH2_100','CCOH-H2O-ele_100','OCCH2_100']:
    species_definitions[sp] = species_definitions['OCCO_100'].copy()

sigma_input = ['CH', 1]
Upzc = 0.00
species_definitions['CO_100']['sigma_params']=[-0.007868567267026672, -0.011350225022150495]
species_definitions['CO-CO_100']['sigma_params']=[0.3473016059671863, -1.65801165572945408]
species_definitions['OCCO_100']['sigma_params']=[0.6138919434717214, -1.7139595521279358]
species_definitions['OCCOH_100']['sigma_params']=[0.48645051933907024, -0.9049509721277826]
species_definitions['HOCCOH_100']['sigma_params']=[0.08611589437959201, 1.0389179051567232]
species_definitions['CCO_100']['sigma_params']=[0.17928131249200074, -0.9859665322047373]
species_definitions['CCOH_100']['sigma_params']=[0.05358796004304156, 1.0505458528246443]
species_definitions['HCCO_100']['sigma_params']=[-0.012044445844595849, 0.22662841070650575]
species_definitions['CCHO_100']['sigma_params']=[0.4433796722495982, -1.3551102625844105]
species_definitions['CC_100']['sigma_params']=[0.26707983721548845, -0.3563150359657281]
species_definitions['H2CCO_100']['sigma_params']=[0.056526361422395116, -0.1298631837225379]
species_definitions['HCCOH_100']['sigma_params']=[0.034029645491422886, -1.1297143663238989]
species_definitions['HCCHO_100']['sigma_params']=[-0.005378573227952107, 0.017474447879992416]
species_definitions['CCH_100']['sigma_params']=[-0.14011016205560217, 0.9122567001114488]
species_definitions['OHCCH2_100']['sigma_params']=[0.23015974211199086, -3.0381698237794756]
species_definitions['H3CCO_100']['sigma_params']=[-0.04213684353513637, 0.0012641742580600246]
species_definitions['H2CCOH_100']['sigma_params']=[0.08774437652694242, 0.34984637183821166]
species_definitions['HCCH_100']['sigma_params']=[-0.11876190528096058, 0.7730448512517086]
species_definitions['OCHCH3_100']['sigma_params']=[-0.1573172511872385, 0.2539067725374835]
species_definitions['H2CCH2O_100']['sigma_params']=[-0.10525316303091975, 0.7301485942106798]
species_definitions['O_100']['sigma_params']=[0.04128867324465629, 0.2367437732076605]
species_definitions['H2CCH_100']['sigma_params']=[-0.16050983422450804, 0.663782643652764]
species_definitions['H3CCH2O_100']['sigma_params']=[-0.27642281719017586, 0.8031121026359245]
species_definitions['OH_100']['sigma_params']=[0.25413030588249225, -0.8363208446700268]

data_file = 'mkm.pkl'
