# The Woolz Image Processing System

## Contents

[Brief Description and History](#Brief-Description-and-History)

[Obtaining Woolz](#Obtaining-Woolz)

[Building Woolz](#Building-Woolz)

[Using Woolz](#Using-Woolz)

[File Formats](#File-Formats)

[Documentation](#Documentation)

[Writing Woolz Software](#Writing-Woolz-Software)


## Brief Description and History

Woolz is a set of software libraries and executables for image processing
and pattern recognition that was initially developed at the MRC Clinical
Population and Cytogenetics Unit. This became the MRC Human Genetics
Unit (HGU) and in 2011 became part of the University of Edinburgh. Woolz
was initially developed for fast microscope slide scanning, chromosome
image analysis, pattern recognition and a wide range of image processing
and analysis problems. The original authors of the software are
Dr. Denis Rutovitz and Dr Jim Piper although
the software has developed and expanded considerably since their initial
work.

Woolz was adopted as the standard for the Mouse Atlas Databases
(http://www.emouseatlas.org) and was used for all the reconstructions,
anatomical, gene-expression and spatial domains. This was because the
interval coding used in Woolz provided significant computing advantages
in a range of image processing functions specifically in set operations
(such as union, intersection, etc.) morphological operations (erosion,
dilation etc.) and other binary image processing such as distance
transforms, segmentation and labelling.
In general Woolz is very efficient with respect to both memory use and time.

The Woolz data structures are, in general, compact and in terms of
grey-level data minimise memory usage without compression (which in
principle could also be applied). It is especially efficient for
morphological and set operations because of the way it's 2 and 3D
domain objects are encoded using intervals.

All pixel or voxel based image processing systems which allow image
query via floating point coordinates need to define what a pixel
is by the relationship between the integer pixel/voxel coordinates and
the floating point query coordinates. In Woolz a pixel or voxel is
defined to be an image sample centred at the integer coordinate.
Following from this definition the C macro WLZ_NINT(), which returns
the nearest integer to the given float, should be used in queries of
pixel/voxel values when using floating point coordinates.

During it's use by the Mouse Atlas Project were many
significant developments in Woolz, primarily focused on
3D reconstruction, transforms, registration and warping.

Since 2019 the development of Woolz has primarily been driven by the
<a name="GCA">
<a
href="https://www.ed.ac.uk/comparative-pathology/the-gut-cell-atlas-project">
Gut Cell Atlas Project</a>
funded by the
<a
href="https://helmsleytrust.org/">
Helmsley Charitable Trust</a>
with significant developments including splines (for midline representation),
JSON encoding and greater use from Python.

Woolz has been written in ANSI standard C so that it will build and
run on all computing platforms that support some basic requirements,
such as supporting at least 32 bit integers and IEEE floating point.
The software is know to build on GNU/Linux systems,
but will also build on MacOS or using MingW on Windows systems.


The code is partitioned into the following modules:

| Module | Description |
| --- | --- |
| libAlc | Library providing generic data structures and memory allocation functions |
| libAlg | Library providing basic numerical algorithms |
| libbibfile  | Library with bibfile style input/output functions |
| libhguDlpList|  Library with a generic doubly linked pointer list |
| libReconstruct| Library with code for 3D alignment of 2D section images to form a 3D image |
| libWlz|         The Woolz image processing library |
| libWlzBnd|      Library with small functions that bind Woolz to other languages |
| libWlzExtFF|    Library for external data format input/output |
| binWlz|         Small command line based Woolz programs |
| binWlzApp|      More small command line based Woolz programs |
| binWlzExtFF|    Small command line based Woolz programs which use external file formats |
| binWlzTst|      Small command line based test programs for Woolz |

The authors include (in sort order):

| <!-- --> |
| --- |
| Bill Hill |
| Christophe Dubreuil |
| Elizabeth Guest |
| Jianguo Rao |
| Jim Piper |
| Konstantinos Liakos |
| Margaret Stark |
| Nick Burton |
| Richard Baldock |

To contact the authors please raise a Github issue.

## Obtaining Woolz

Woolz is available as source from GitHub (https://github.com/ma-tech/)
and can be downloaded pre-built for some architectures.

## Building Woolz

Woolz should build easily on most modern systems that have the GNU Build
system (see http://en.wikipedia.org/wiki/GNU_build_system).

In most cases the simple script build.sh should be sufficient to build Woolz.
This script simply runs autoreconf followed by ./configure with various
options. It's probably best to copy build.sh to mybuild.sh and then set the
required options for the build you want.

Sometimes the configuration will fail with messages relating to libtool.
In these cases running:
  autoreconf -i --force
may fix the problems.

Complains from configure that m4 can't be found can probably be fixed
by running:
  automake --add-missing

The build script build.sh simply runs autoreconf followed by ./configure
with the appropraite options.

A prefix can be given in the configure stage to define where the programs,
libraries are will be installed, conventionally this is /opt/MouseAtlas.
Use

./configure --help

to see all the options available.

## Using Woolz

There are over 200 small command line programs within Woolz. All
of these accept -h as an argument to show their usage. On unix-like
systems it is common to combine these small programs into a single
command line with pipes. These programs include those for:
  applying affine, basis function and mesh based transforms
  registration of spatial domain objects and surfaces
  reconstruction from serial sections
  morphological operations (erosion, dilation, etc) and fill
  set operations (union, intersection etc)
  feature extraction
  mesh generation
  contour generation
  convex hulls
  distance transforms
  histograms
  thresholding and labelling
  grey and colour image value filters

As an example, the following thresholds a 3D image, applies erosion and
dilation to remove small isolated regions, labels (segments) the image
into separate objects and the prints the volume of each isolated object:

```
prompt% WlzThreshold -v135 -L ts14.wlz | \
        WlzErosion  -c26 -r2 | \
        WlzDilation -c26 -r2 | \
	WlzLabel | \
	WlzVolume
Object 1: number of voxels = 66
Object 2: number of voxels = 45
Object 3: number of voxels = 45
Object 4: number of voxels = 66
Object 5: number of voxels = 159
Object 6: number of voxels = 45
Object 7: number of voxels = 45
Object 8: number of voxels = 58
Object 9: number of voxels = 276
Object 10: number of voxels = 45
Object 11: number of voxels = 138
Object 12: number of voxels = 45
Object 13: number of voxels = 45
Object 14: number of voxels = 137
Object 15: number of voxels = 45
```

Woolz can also be used as a set of C libraries, or via a binding
to another language. Currently these bindings exist for Java, Python
and R. For details of these bindings see the JavaWoolz, PyWoolz and
RWoolz repositories.

## File FormatS

The encoding of Woolz objects when serialised to files is defined only
by the source code in the WlzReadObj()/WlzWriteObj() functions. For
historical reasons there is no unique identifier (magic number) for
Woolz objects. Woolz objects can also (usualy with some loss of information)
be written to other (external) file formats.

The program WlzExtFFConvert can be used to convert between supported
file formats. As with all the small Woolz programs, the -h option will
show usage, but for WlzExtFFConvert it will also list the file formats
which are understood:

```
prompt% WlzExtFFConvert -h
Usage: WlzExtFFConvert [-h] [-s] [-b<background>]
                       [-d<min-dimension>] [-D<max-dimension>]
                       [-f<input format>] [-F<output format>]
                       [-x<x size>] [-y<y size>] [-z<z size>]
                       [-o<output file>] [<input file>)]
Converts objects between one file format and another, neither of
which need be the Woolz data file format.
Version: 1.8.2
Options:
  -h    Help, prints this usage information.
  -s    Split labeled volumes into domains.This will also split
        a tiled or pyramidal tiff into resolution slices and tiles.
  -G    Apply grey value transforms.
  -S    Use spatial transforms to create WLZ_TRANS_OBJ objects
        (by default just offsets are applied).
  -b#   Set background to value,
  -d#   Set size of minimum dimension, i.e. min of width or height
  -D#   Set size of maximum dimension, i.e. max of width or height
  -f#   Input file format.
  -F#   Ouput file format.
  -o#   Output file name.
  -x#   X voxel/pixel size.
  -y#   Y voxel/pixel size.
  -z#   Z voxel/pixel size.
The known file formats are:
  Description                                       Extension
  ***********                                       *********
  Amira Lattice                                     am
  Microsoft Bitmap                                  bmp
  Stanford Density                                  den
  Netgen neutral mesh format                        emt
  Graphics Interchange Format                       gif
  ANALYZE HDR                                       hdr
  Image Cytometry Standard                          ics
  IPLab                                             ipl
  JPEG                                              jpg
  JSON encoded Woolz                                jsn
  Pascal Frey's medit tetrahedral mesh format       mesh
  Neuroimaging Informatics Technology Initiative    nii
  Jonathan Shewchuk's mesh format                   node
  Utah nearly raw raster data (NRRD) format         nrrd
  Wavefront                                         obj
  BioRad Confocal                                   pic
  Riken PLY2                                        ply2
  PNM                                               pnm
  Drishti dot NC format                             pvl.nc
  Raw                                               raw
  SLC                                               slc
  GRUMMP SMESH                                      smesh
  Stereolithography format                          stl
  Tiff                                              tif
  Text                                              txt
  Sunvision VFF                                     vff
  GRUMMP VMESH                                      vmesh
  Visualization Toolkit VTK                         vtk
Simple example:
  WlzExtFFConvert -f wlz -F slc <in.wlz >out.slc
  Converts the Woolz object in.wlz to an SLC data file out.slc
More complex example:
  WlzExtFFConvert -f den -F pnm -o out.pgm in.den
  Converts the Stanford density file in.den to a series of PGM files
  each with a name of the form out000001.pgm where the number
  encodes the image plane and a control file which specifies the
  volume origin, size and voxel dimensions.
By default objects are read from the standard input and written to
the standard output.
File formats which use more than one file can not be read or written
using the standard input or standard output.
The TIFF file format must be read/written from/to a file i.e. not
from/to stdin or stdout
Not all formats can retain position information i.e. they can
only keep the size of the bounding box. In these formats the
size of the bounding box is maintained but the position is set to
(0,0). This implies that conversion back to woolz will, in general,
result in a shifted image, i.e. registration is lost. Most 3D
formats encode this data, of the 2D formats only woolz can retain
all offsets, TIFF can only encode positive offsets.
```

The library functions WlzEffReadObj() and WlzEffWriteObj() (in libWlzExtFF)
can also be used to read and write non-Woolz format files.

##  Documentation

Woolz uses Doxygen (http://www.doxygen.org) for documentation,
although the best documentation is (as always) the source code itself,
the Doxygen documentation is available from
https://ma-tech.github.io/Woolz/documentation/html_Core/index.html .

## Writing Woolz Software

The best way to understand Woolz is through the source code, but what
follows attempts to give an overview and may be some help.

### The Woolz Object

```c
  typedef struct _WlzObject
  {
    WlzObjectType      type;
    int                linkcount;
    WlzDomain          domain;
    WlzValues          values;
    WlzPropertyList    *plist;
    struct _WlzObject  *assoc;
  } WlzObject;
```

The fields encode the type, link count, spatial domain,
values, properties and any associated objects of an object.
The type simply encodes what the object is (3D image, 2D polygon, ...).
Object use reference counting via the linkcount makes many Woolz
operations, such as thresholding, extremely efficient with regard to
both space and time. The link (reference) count of an object is
incremented when an object is assigned using WlzAssignObject() and
decremented when an object is freed using WlzFreeObj(). Typical
usage looks like:

```c
  WlzObject *obj;

  obj = WlzAssignObject(WlzReadObj(stdin, &errNum), NULL);
  /* Do something with obj */
  (void )WlzFreeObj(obj);
```

The domain of an object is the spatial extent within which the object
is defined, some other spatial representation (including meshes), some
transformation or bizarrely a histogram.

The values of an object are some values, such as intensities, which are
only defined within the objects domain.

There is a core Woolz object type:

```c
  typedef struct _WlzCoreObject
  {
    WlzObjectType type;
    int           linkcount;
  } WlzCoreObject;
```

which is sufficient to determine the type of and object and to either
assign it or free it. So, "there's inheritance too, sort of". The
other type of top level Woolz object is an array of objects:

```c
  typedef struct _WlzCompoundArray
  {
    WlzObjectType type;
    int           linkcount;
    WlzObjectType otype;
    int           n;
    WlzObject     **o;
    WlzPropertyList *plist;
    WlzObject     *assoc;
  } WlzCompoundArray;
```

This is similar to the main Woolz object type but has an array of objects
rather than a domain and values (variants of this compound object exist
including variants with linked lists).

#### Woolz Domains

The domain of an object is (in most cases) the spatial description of
and object. For 2D or 3D domain objects, such as the EMAGE anatomy
domains or the EMAGE model embryos, the domain encodes the region of
space that the anatomical or embryo component occupies.
Just as there is a core object there is a core domain which has the
same uses as the core object:

```c
  typedef struct< _WlzCoreDomain
  {
    WlzObjectType   type;
    int             linkcount;
    void            *freeptr;
  } WlzCoreDomain;
```

The additional field (compared to WlzCoreObject) is the free pointer.
This is a pointer to a stack of memory blocks that have been allocated
for the domain. The Woolz Domain is a union of the possible domains:

```c
  typedef union _WlzDomain
  {
    struct _WlzCoreDomain      *core;
    struct _WlzIntervalDomain  *i;
    struct _WlzPlaneDomain     *p;
    struct _WlzPolygonDomain   *poly;
    struct _WlzBoundList       *b;
    struct _WlzHistogramDomain *hist;
    struct _WlzRect            *r;
    struct _WlzFRect           *fr;
    struct _WlzAffineTransform *t;
    struct _WlzWarpTrans       *wt;
    struct _WlzContour         *ctr;
    struct _WlzMeshTransform   *mt;
    struct _WlzLBTDomain2D     *l2;
    struct _WlzLBTDomain3D     *l3;
    struct _WlzCMesh2D         *cm2;
    struct _WlzCMesh2D5        *cm2d5;
    struct _WlzCMesh3D         *cm3;
    struct _WlzPoints          *pts;
    struct _WlzLUTDomain       *lut;
    struct _WlzThreeDViewStruct *vs3d;
  } WlzDomain;
```

These include domains for 2D and 3D spatial regions (WlzIntervalDomain
and WlzPlaneDomain) as well as polygons, boundaries, contours, meshes
and transforms (such as affine, basis function and mesh transforms).
Transforms are domains since they are spatial mappings. Histograms
are domains too (for historical reasons).

#### Woolz Values

The values of an object (again in most cases) represent the actual
values that are embedded in the space defined by the domain. An
object with a domain but without values is perfectly valid and are
frequently used that way, eg to represent some anatomical region.
The core values type is:

```c
  typedef struct _WlzCoreValues
  {
    WlzObjectType type;
    int       linkcount;
  } WlzCoreValues;
```

The Woolz Values is (similarly to the Woolz Domain) a union of the
possible Woolz values:

```c
  typedef union _WlzValues
  {
    struct _WlzCoreValues     *core;
    struct _WlzRagRValues     *v;
    struct _WlzRectValues     *r;
    struct _WlzIntervalValues *i;
    struct _WlzConvHullValues *c;
    struct _WlzVoxelValues    *vox;
    struct _WlzObject         *obj;
    struct _WlzFeatValues     *fv;
    struct _WlzRectFeatValues *rfv;
    struct _WlzIndexedValues  *x;
    struct _WlzTiledValues    *t;
    struct _WlzLUTValues      *lut;
  } WlzValues;
```

The types of values include image values (WlzRagRValues, WlzRectValues,
WlzIntervalValues and WlzVoxelValues), convex hulls, indexed values
(used with meshes) look up tables and features. The top level Woolz object
itself is also a valid values! This allows objects to be defined which
include a spatial mapping, such as an affine transform, without computing
the transformation of the object.

#### Woolz Properties

Property lists are arbitrary lists of properties that are associated with
an object such as a meaningful name. They were once used extensively,
but have now been reinstated as genuine linked lists. Properties may
be used to record the history of an object or object names.

### Access Methods

Given that Woolz objects are complex data structures, access functions
are required to use their elements. These access methods are for spatial
domain objects with values (images) but simple adaptations can make the
code applicable to spatial domain objects without values.

#### Interval Scanning

Processing 2D Woolz domain objects is most efficiently done by scanning
through objects using blocks of contiguous pixels.

```c
  WlzGreyP gP,
  WlzObject *obj;
  WlzIntervalWSpace iWsp;
  WlzGreyWSpace gWsp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  errNum = WlzInitGreyScan(obj, &iWsp, &gWsp);
  while((errNum == WLZ_ERR_NONE) &&
        ((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE))
  {
    gP = gWsp.u_grintptr;
    switch(gWsp.pixeltype)
    {
      case WLZ_GREY_INT:
        for(iPos = iWsp.lftpos; iPos &lt;= iWsp.rgtpos; ++iPos)
        {
          *(gP.inp + iPos) /= 2;
        }
        break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
```

#### Random Access

There are random access functions for establishing whether some vertex
within a spatial domain object (WlzInsideDomain()) and the value at
some position in space, either inside or outside the object:

```c
  int       val;
  WlzIVertex3 pos;
  WlzObject *obj;
  WlzGreyValueWSpace *gVWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    WlzGreyValueGet(gVWSp, pos.vtZ, pos.vtY pos.vtX);
    switch(gVWSp-&gt;gType)
    {
      case WLZ_GREY_INT:
        val = (*(gVWSp->gVal)).inv;
        break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
        break;
    }
  }
```

#### Iterators

Simple iterators are another way to access the values of spatial domain
objects in scan order. Unlike WlzNextGreyInterval() the iterators work
for both 2D and 3D. In the example below all (integer) values of the
spatial domain object are incremented from 0.

```c
  WlzErrorNum SetInvGreyValues(WlzObject *obj)
  {
    int         i = 0;
    WlzIterateWSpace *itWSp = NULL;
    WlzErrorNum   errNum = WLZ_ERR_NONE;

    itWSp = WlzIterateInit(obj, WLZ_RASTERDIR_IPILIC, 1, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      while((errNum = WlzIterate(itWSp)) == WLZ_ERR_NONE)
      {
	*(itWSp->gP.inp) = i++;
      }
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
    }
    WlzIterateWSpFree(itWSp);
    return(errNum);
  } 
```


### Geometric Models

Woolz geometric models provide a unified representation of both 2D and
3D geometric models composed of simplices. They are capable of simple planar
straight line graphs in 2D and both manifold and non-manifold surfaces
(with non intersecting elements) in 3D.

### Contours

Woolz geometric models are used to represent contours such as iso-value
contours extracted from domain objects.

### Affine Transforms

An affine transform is a transformation which preserves lines and the
parallelism of lines, but not necessarily lengths or the angles between
(non parallel) lines. Affine transforms in Woolz are stored as homogeneous
3x3 and 4x4 matrices (actually 4x4 arrays are used for both but 2D
transforms access them as if they are 3x3).

### Basis Function Transforms

Basis function transforms allow displacements at discrete points to be
interpolated throughout an object's domain. These are the sum of radially
symmetric component transforms, with each component transform centred on
one of the points.

### Meshes

Data structures and methods exist within Woolz for both convex and
conforming meshes composed of simplices in 2, 2.5 (2D topology but
3D geometry, ie surfaces in 3D space) and 3D. Transforms may be
built using the meshes by associating  displacements with the mesh
nodes, in the case of conforming meshes this is done using indexed
values.

