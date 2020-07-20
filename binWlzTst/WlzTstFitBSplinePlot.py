###
#ident "University of Edinburgh $Id$"
#
# \file         WlzTstFitBSplinePlot.py
# \author       Bill Hill
# \date         June 2020
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2020],
# The University Court of the University of Edinburgh,
# Old College, Edinburgh, UK.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief    Simple plotting script for evaluation output from
#           WlzTstFitBSpline
# \ingroup  BinWlzTst
###

import sys
import json
import math as m
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) == 2:
  jf = sys.argv[1]
else:
  jf = 'out.jsn'

with open(jf, 'r') as f:
  t=f.read()

ntsz = 0.1
test = json.loads(t)
dim = test["outdata"]["dim"]
sm = test["outdata"]["smooth"]
k = test["outdata"]["order"]
n = test["outdata"]["knots"]
din = np.array(test["indata"]["points"])
dot = np.array(test["outdata"]["points"])
try:
  normalised = test["outdata"]["normalised"]
except:
  normalised = False
#end
try:
  postnt = test["postnt"]
except:
  postnt = None
#end

fig = plt.figure()
if(dim == 2):
  ax = plt.subplot(111)
  ax.plot(din.transpose()[0], din.transpose()[1],
          'k--', label='Input data', marker='.', markerfacecolor='k')
  ax.plot(dot.transpose()[0], dot.transpose()[1],
          'b', linewidth=2.0, label='Interpolated B-spline')
  if(normalised):
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
  #end
  if(bool(postnt)):
    nln = np.zeros((2, 2))
    pnt = np.array([postnt[0], postnt[1]])
    tnt = np.array([postnt[2], postnt[3]])
    if(m.fabs(tnt[0]) > m.fabs(tnt[1])):
      nrm = np.array([-tnt[1] / tnt[0], 1.0])
    else:
      nrm = np.array([1.0, -tnt[0] / tnt[1]])
    #end
    nrm = nrm / np.linalg.norm(nrm)
    tgt = np.zeros((2, 2))
    for i in range(0, 2):
      tgt[i][0] = pnt[i]
      tgt[i][1] = pnt[i] + ntsz * tnt[i]
    #end
    for i in range(0, 2):
      nln[i][0] = pnt[i] - ntsz * nrm[i]
      nln[i][1] = pnt[i] + ntsz * nrm[i]
    #end
    ax.plot(tgt[0], tgt[1], 'magenta', linewidth=2.0, label='Tangent')
    ax.plot(nln[0], nln[1], 'red', linewidth=2.0, label='Normal')
  #end
else:
  ax3d = fig.add_subplot(111, projection='3d')
  #  ax3d.plot(din.transpose()[0], din.transpose()[1], din.transpose()[2],
  #          'k--', label='Input data', marker='.', markerfacecolor='k')
  ax3d.plot(dot.transpose()[0], dot.transpose()[1], dot.transpose()[2],
            'b', linewidth=2.0, label='Interpolated B-spline')
  if(normalised):
    ax3d.set_xlim(-0.1, 1.1)
    ax3d.set_ylim(-0.1, 1.1)
    ax3d.set_zlim(-0.1, 1.1)
  #end
  if(bool(postnt)):
    nln = np.zeros((3, 5))
    pnt = np.array([postnt[0], postnt[1], postnt[2]])
    tnt = np.array([postnt[3], postnt[4], postnt[5]])
    # Find max component of tangent
    mi = 0;
    for i in range(1, 3):
      if(m.fabs(tnt[i]) > m.fabs(tnt[mi])):
        mi = i
      #end
    #end
    # Find vectors orthogonal to the tangent
    t0 = np.array([tnt[0], tnt[1], tnt[2]])
    t0[mi] = 0.0
    t0 = t0 / np.linalg.norm(t0)
    t1 = np.cross(tnt, t0)
    t1 = t1 / np.linalg.norm(t1)
    t2 = np.cross(tnt, t1)
    t2 = t2 / np.linalg.norm(t2)
    for i in range(0, 3):
      nln[i][0] = pnt[i] - ntsz * t1[i] - ntsz * t2[i]
      nln[i][1] = pnt[i] + ntsz * t1[i] - ntsz * t2[i]
      nln[i][2] = pnt[i] + ntsz * t1[i] + ntsz * t2[i]
      nln[i][3] = pnt[i] - ntsz * t1[i] + ntsz * t2[i]
      nln[i][4] = nln[i][0]
    #end
    t0 = np.zeros((3, 2))
    for i in range(0, 3):
      t0[i][0] = pnt[i]
      t0[i][1] = pnt[i] + ntsz * tnt[i]
    #end
    tgt = t0
    ax3d.plot(tgt[0], tgt[1], tgt[2], 'magenta', linewidth=2.0, label='Orange')
    ax3d.plot(nln[0], nln[1], nln[2], 'red', linewidth=2.0, label='Normal')
  #end
#end
#plt.legend(['Input data', 'Interpolated B-spline', 'Tangent', 'Normal'],
#           loc='upper center', bbox_to_anchor=(0.05, 0.05))
plt.title('B-Spline interpolation (s=' + 
    str(sm) + ', k=' + str(k) + ', n=' + str(n) + ')')
plt.savefig('plot.png')
#plt.show()

