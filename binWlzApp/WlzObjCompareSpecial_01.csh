#! /bin/tcsh -f

set filelist = $1

echo " "
echo Calculating domain comparison values for $filelist
echo Path: `pwd`
echo User: `whoami`
echo Date: `date`
echo " "
echo Output from script: $0

set files = `cat $1`
set i = 1

echo " "
echo Domain 1, Domain 2, volume 1, volume 2, mass 1, mass 2, cm-cm, cm-surf, surf-cm, surf-surf

while ( $i <  $#files )
    set cosmidfile = $files[$i]
    @ i += 1
    set terrfile = $files[$i]
    @ i += 1
    echo -n $cosmidfile, $terrfile, " "
    cat $cosmidfile $terrfile | WlzObjCompareSpecial_01
end
