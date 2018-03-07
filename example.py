from ssf import rdf, structureFactor, structureFactorPowder
from ssf.RDF.RDF import ssf

import numpy as np
import itertools as it
import pylab as pl
import time

def generate_sc(size, n):
  """
  Generate the positions of a simple cubic crystal in a box of
  length size with n atoms in each direction (order parameter =
  size/n)
  """

  natoms = n**3
  pos = range(n)
  x = np.zeros((natoms, 3))
  i = 0
  for px, py, pz in it.product(pos, pos, pos):
    x[i] = (px, py, pz)
    x[i] = x[i] * (size/n)
    i += 1
  return x

#pl.matplotlib.interactive(True)
ldorder = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,
           350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,
           3074, 3890, 4334, 4802, 5294, 5810]
#ldorder = [5810]
lpart = range(2, 10)
time_powder = []
slope_leb = []
for i in lpart:
  size = float(i)
  k = np.linspace(0, 20.0, 100)
  x = generate_sc(1.0, i)
  natoms = np.shape(x)[0]
  b = time.time()
  exact = structureFactorPowder(x, size, k)
  t_pow = time.time() - b
  print "Elapsed in ssf powder: {0}s".format(t_pow)
  time_powder.append(t_pow)
  fig, ax = pl.subplots()


  mean = []
  maximum = []
  t_leb = []
  leb = [6, 50, 194, 590, 3890, 5810]
  for l in ldorder:
    b = time.time()
    s3 = structureFactor(x, size, k, rep=1, lebedev=l)
    td = time.time() - b
    print "Elapsed in ssf leb({0}): {1}s".format(l, td)
    t_leb.append(td)
    mean.append(np.mean(abs(exact[:, 1]-s3[:, 1])))
    maximum.append(np.max(abs(exact[:, 1]-s3[:, 1])))
    if l in leb:
      ax.plot(k, s3[:, 1], label='Leb - {0}'.format(l))
  ax.plot(k, exact[:, 1], label='Exact')
  ax.set_ylabel('S(q)')
  ax.set_xlabel('q')
  ax.legend()
  ax.set_title('Sq para {0}'.format(i))
  fig.tight_layout()
  fig.savefig('{0}_sq.png'.format(i))
  pl.close()
  
  fig, ax2 = pl.subplots()
  ax2.plot(ldorder, mean, '-o', label='Mean')
  ax2.plot(ldorder, maximum, '-o', label='Max')
  ax2.set_xscale("log")
  ax2.set_yscale("log")
  ax2.set_ylabel('Error absoluto')
  ax2.set_xlabel('Numero de puntos')
  ax2.legend()
  ax2.set_title('Error para {0}'.format(i))
  fig.tight_layout()
  fig.savefig('{0}_err.png'.format(i))
  pl.close()
  
  fig, ax3 = pl.subplots()
  ax3.plot(ldorder, t_leb, '-o')
  ax3.set_ylabel('Tiempo')
  ax3.set_xlabel('Numero de puntos')
  pol = np.polyfit(ldorder, t_leb, 1)
  ax3.plot(ldorder, np.polyval(pol, ldorder))
  ax.set_title('Tiempo para {0} ({1})'.format(i, pol[0]))
  fig.tight_layout()
  fig.savefig('{0}_time.png'.format(i))
  pl.close()

  slope_leb.append(pol[0])

