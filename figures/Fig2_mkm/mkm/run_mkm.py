from catmap import ReactionModel
from catmap import analyze
from string import Template
import os,sys
import pickle

model = ReactionModel(setup_file = 'reactions.mkm')
#model.output_variables+=['production_rate', 'free_energy','rate_control','selectivity_control']
model.run()

exit()
vm = analyze.VectorMap(model)
ma = analyze.MechanismAnalysis(model)
ma.energy_type = 'free_energy'
label_size = 10

vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e+10
fig = vm.plot(save=False)
fig.savefig('output/rate.pdf')

vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.colorbar = True
vm.min = 1e-5
vm.max = 1e+7
fig = vm.plot(save=False)
fig.savefig('output/production_rate.pdf')

vm = analyze.VectorMap(model)
vm.log_scale = False
vm.unique_only = False
vm.plot_variable = 'coverage'
vm.min = 0
vm.max = 1
fig = vm.plot(save=False)
fig.savefig('output/coverage.pdf')

vm = analyze.VectorMap(model)
vm.log_scale = True
vm.unique_only = False
vm.plot_variable = 'coverage'
vm.min = 1e-20
vm.max = 1
fig = vm.plot(save=False)
fig.savefig('output/coverageLog.pdf')

mm = analyze.MatrixMap(model)
mm.plot_variable = 'rate_control'
mm.log_scale = False
mm.min = -2
mm.max = 2
mm.plot(save='output/rate_control.pdf')

mm = analyze.MatrixMap(model)
mm.plot_variable = 'selectivity_control'
mm.log_scale = False
mm.min = -2
mm.max = 2
mm.plot(save='output/selectivity_control.pdf')


