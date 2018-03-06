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

def myssf(x, size, q, pbc=False):
  """
  From a series of positions x in a cubic box of length size we get
  the structure factor for momentum q
  """

  natoms = np.shape(x)[0]
  sf = 0.0
  for i in range(natoms):
    x1 = x[i]
    for j in range(i+1, natoms):
      x2 = x[j]
      dx = x2 - x1
      if pbc:
        for i in range(3):
          if dx[i] >  size/2: dx[i] -= size
          if dx[i] < -size/2: dx[i] += size
      r = np.linalg.norm(dx)
      sf += 2*np.sin(q*r)/(q*r)
      sf /= natoms
  sf += 1
  sf[q==0] = natoms
  return sf

size = 1.0
x = generate_sc(size, 10)
natoms = np.shape(x)[0]

"""r = rdf(x, size, 1000, pbc=False)
s = ssf(r, natoms/size**3, pbc=False)
r = rdf(x, size, 1000, pbc=True)
s2 = ssf(r, natoms/size**3, pbc=True)
k = s[:, 0]
k = k[k < 20.0]
fig, ax = pl.subplots()
ax.plot(s2[:, 0], s2[:, 1], label='PBC')
ax.plot(s[:, 0], s[:, 1], label='No PBC')
"""
k = np.linspace(0, 20.0, 100)
print "Exact - C"
b = time.time()
exact = structureFactorPowder(x, size, k)
t_pow = time.time() - b
print "Elapsed in ssf powder: {0}s".format(t_pow)

#print "Exact"
#exact = myssf(x, size, k)
#fig, ax = pl.subplots()
#ax.plot(k, exact, label='Exact')
#ax.plot(k, exact_c[:, 1], label='Exact in C')
#pl.show()

ldorder = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,
           350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,
           3074, 3890, 4334, 4802, 5294, 5810]
fig, ax = pl.subplots()


d = []
t_leb = []
leb = [6, 50, 194, 590, 3890, 5810]
for l in ldorder:
  print "Calculating {0}".format(l)
  b = time.time()
  s3 = structureFactor(x, size, k, rep=1, lebedev=l)
  td = time.time() - b
  print "Elapsed in ssf leb({0}): {1}s".format(l, td)
  t_leb.append(td)
  d.append(np.mean(abs(exact[:, 1]-s3[:, 1])))
  if l in leb:
    ax.plot(k, s3[:, 1], label='Leb - {0}'.format(l))
ax.plot(k, exact[:, 1], label='Exact')
ax.legend()

fig, ax2 = pl.subplots()
ax2.plot(ldorder, d, '-o')
ax2.set_xscale("log")
ax2.set_yscale("log")

fig, ax3 = pl.subplots()
ax3.plot(ldorder, t_leb, '-o')
