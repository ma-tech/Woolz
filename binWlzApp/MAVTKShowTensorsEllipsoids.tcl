#!/usr/bin/wish
##
# \file         MAVTKShowTensorsEllipsoids.tcl
# \author       Bill Hill
# \date         February 2013
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2013],
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
# \ingroup      BinWlzApp
# \brief        Display tensor ellipsiods using VTK.
##

package require vtk
package require vtkinteraction

source /opt/vis/share/VTK/colors.tcl

set inFileName "out.vtk"

vtkPolyDataReader reader
  reader SetFileName $inFileName
  reader SetTensorsName "tensors"
  reader Update

vtkSphereSource sphere
  sphere SetRadius 16.0
  sphere SetThetaResolution 12
  sphere SetPhiResolution 12

vtkAxes axes
  axes SetScaleFactor 16.0

vtkTubeFilter tubeAxes
  tubeAxes SetInput [axes GetOutput]
  tubeAxes SetRadius 2.0
  tubeAxes SetNumberOfSides 4

vtkTensorGlyph ellipsoids
  ellipsoids SetInput [reader GetOutput]
  ellipsoids SetSource [sphere GetOutput]
#  ellipsoids SetSource [tubeAxes GetOutput]
  ellipsoids SetScaleFactor 0.2
#  ellipsoids ClampScalingOn
#  ellipsoids ColorGlyphsOff
#  ellipsoids ColorGlyphsOn
#  ellipsoids ExtractEigenvaluesOn
#  ellipsoids SetColorModeToEigenvalues
#  ellipsoids ThreeGlyphsOff
   ellipsoids SetColorModeToEigenvalues
#   ellipsoids ThreeGlyphsOn
   ellipsoids SetMaxScaleFactor 8.0

vtkPolyDataNormals ellipNormals
  ellipNormals SetInput [ellipsoids GetOutput]

vtkPolyDataMapper ellipMapper
    ellipMapper SetInput [ellipNormals GetOutput]

vtkActor ellipActor
    ellipActor SetMapper ellipMapper

vtkRenderer ren

vtkRenderWindow renWin
  renWin AddRenderer ren
  renWin SetSize 512 512
  renWin StereoCapableWindowOn
  renWin SetStereoTypeToCrystalEyes

vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin
  iren AddObserver UserEvent {wm deiconify .vtkInteract}
  iren AddObserver ExitEvent {exit}

  ren AddActor ellipActor
  ren SetBackground 255 255 255

  renWin Render

wm withdraw .

tkwait window .
