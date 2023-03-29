import pylab as pl
import matplotlib
import numpy as np
from matplotlib.colors import LogNorm

font = {'family' : 'Times New Roman',
        'weight' : 'regular',
        'size'   : 35}

matplotlib.rc('font', **font)


a = np.array([[1e-5,1e5]])
pl.figure(figsize=(9, 1))
img = pl.imshow(a, cmap="jet",norm=LogNorm())
pl.gca().set_visible(False)
cax = pl.axes([0.1, 0.7, 0.8, 0.2])
pl.colorbar(orientation="horizontal", cax=cax,norm=LogNorm())
pl.savefig("colorbar_current.pdf")
#pl.tight_layout()
pl.show()
pl.close()

a = np.array([[0,1]])
pl.figure(figsize=(9, 1))
img = pl.imshow(a, cmap="RdYlGn")
pl.gca().set_visible(False)
cax = pl.axes([0.1, 0.7, 0.8, 0.2])
pl.colorbar(orientation="horizontal", cax=cax)
pl.savefig("colorbar_selectivity.pdf")
pl.tight_layout()
pl.show()
pl.close()


