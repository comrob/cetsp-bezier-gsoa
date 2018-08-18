#!/bin/sh

if [ `uname` = FreeBSD ]
then
   cmd=gmake 
else
   cmd=make 
fi

#export CXXFLAGS=-std=c++17

rm -rf build

mkdir -p build


cd build && cmake -D CMAKE_INSTALL_PREFIX=. ..
$cmd -C crl install
$cmd cetsp-bezier-gsoa
$cmd -C cetsp-bezier-gsoa install
cd -

