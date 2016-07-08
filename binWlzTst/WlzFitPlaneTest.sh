#!/bin/sh
#set -x

PIT=0
YAW=0
D=0
FX=0
FY=0
FZ=0

if [ $# -ge 1 ]
then
  PIT=$1
  if [ $# -ge 2 ]
  then
    YAW=$2
  fi
fi

PIT=`echo "$PIT * 3.141592653 / 180" | bc -l`
YAW=`echo "$YAW * 3.141592653 / 180" | bc -l`
# echo "--- $PIT --- $YAW ---"

rm -f j0.wlz

for P in '181,62,39' \
         '105,192,97' \
	 '281,185,282' \
	 '342,81,235' \
         '118,139,54' \
	 '221,146,173' \
	 '312,147,274' \
	 '225,69,94' \
         '190,148,142' \
	 '141,208,154'
do
  # echo "--- $P ---"
  X=`echo $P | awk -F',' '{print $1}'`
  Y=`echo $P | awk -F',' '{print $2}'`
  Z=`echo $P | awk -F',' '{print $3}'`
  if [ $# -eq 0 ]
  then
    T="$X $Y $Z"
  else
    T=`WlzMakeRect -2 -x $X,$X -y $Y,$Y | \
       Wlz3DViewTransformObj -a $PIT,$YAW -f $FX,$FY,$FZ -d $D |
       WlzBoundingBox`
  fi
  # echo "--- $T ---"
  X=`echo $T | awk '{print $1}'`
  Y=`echo $T | awk '{print $2}'`
  Z=`echo $T | awk '{print $3}'`
  echo "$X $Y $Z"
  WlzMakeRect -3  -x $X,$X -y $Y,$Y -z $Z,$Z >> j0.wlz
done

WlzUnion j0.wlz | WlzDilation -c26 -r2 | WlzGreySetValue -g 100 >j1.wlz
