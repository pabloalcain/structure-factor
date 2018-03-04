"""
Calculation of the structure factor
"""

import ctypes as ct
import numpy as np
import os

_DIRNAME = os.path.dirname(__file__)
libssf = ct.CDLL(os.path.join(_DIRNAME, 'libssf.so'))
ssf_c = libssf.ssf

def structureFactor(x, box, k, rep=2, lebedev=194):
  """
  Calculate structure factor.

  Parameters
  ----------

  x : numpy float64 array
      Positions of the particles in the system

  box : numpy float64 array
      Box

  k : numpy array
      Wavenumbers to calculate

  repetitions : int, optional
      Number of repetitions of the principal cell to
      consider. Default value is 2

  lebedev : int, optional
      Number of points in the sphere in which to calculate the
      Lebedev quadrature. Default value is 194

  Returns
  -------

  ssf : dict
      An array with the information of the compute. The first column
      is 'r' and the rest is the structure factor calculated for the
      pair list.
  """

  size_x = box[0][1] - box[0][0]
  size_y = box[1][1] - box[1][0]
  size_z = box[2][1] - box[2][0]
  if size_x != size_y or size_y != size_z:
    raise ValueError("The box should be cubic for this to work")
  else:
    size = size_x
  natoms = np.shape(x)[0]
  npoints = len(k)
  tmp = (ct.c_double * (npoints * 2))()
  x_p = x.ctypes.data_as(ct.c_void_p)
  k_p = k.ctypes.data_as(ct.c_void_p)
  ssf_c.argtypes = [ct.c_void_p, ct.c_int, ct.c_double, ct.c_int,
                    ct.c_int, ct.c_int, ct.c_void_p, ct.c_void_p]
  ssf_c(x_p, natoms, size, npoints, lebedev, rep, k_p, tmp)
  ssf = np.frombuffer(tmp, dtype=np.double, count=npoints * 2)
  return ssf.reshape((npoints, 2))
