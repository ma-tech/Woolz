#!/bin/sh

swig -D__GNUC__ -Wall -java -ljava -module Wlz -includeall -ignoremissing -I../libAlc -I../libAlg Wlz.h &&
mkdir -p build && /opt/java/bin/javac -d build *.java &&
gcc -I. -o libWlz_wrap.so -I../libAlc -I../libAlg -fPIC -shared Wlz_wrap.c
