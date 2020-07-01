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

test = json.loads(t)
dim = test["outdata"]["dim"]
sm = test["outdata"]["smooth"]
k = test["outdata"]["order"]
n = test["outdata"]["knots"]
din = np.array(test["indata"]["points"])
dot = np.array(test["outdata"]["points"])

fig = plt.figure()
if(dim == 2):
  ax = plt.subplot(111)
  ax.plot(din.transpose()[0], din.transpose()[1],
          'k--', label='Input data', marker='.',markerfacecolor='k')
  ax.plot(dot.transpose()[0], dot.transpose()[1],
          'b',linewidth=2.0,label='Interpolated B-spline')
else:
  ax3d = fig.add_subplot(111, projection='3d')
  ax3d.plot(din.transpose()[0], din.transpose()[1], din.transpose()[2],
            'k--', label='Input data', marker='.',markerfacecolor='k')
  ax3d.plot(dot.transpose()[0], dot.transpose()[1], dot.transpose()[2],
            'b',linewidth=2.0,label='Interpolated B-spline')
#fi
plt.legend(['Input data', 'Interpolated B-spline', 'True'],
           loc='upper center', bbox_to_anchor=(0.05, 0.05))
plt.title('B-Spline interpolation (s=' + 
    str(sm) + ', k=' + str(k) + ', n=' + str(n) + ')')
plt.show()

