#!/bin/bash

ssmem=https://github.com/LPD-EPFL/ssmem
ssmem_folder=./external/ssmem;
sspfd=https://github.com/trigonak/sspfd
sspfd_folder=./external/sspfd;

if [ ! -d $ssmem_folder ];
then
    git clone $ssmem $ssmem_folder
fi;

if [ ! -d $sspfd_folder ];
then
    git clone $sspfd $sspfd_folder
fi;

mkdir external/lib &> /dev/null;

cd $ssmem_folder;
git pull;
make libssmem.a;
cp libssmem.a ../lib;
cp include/ssmem.h ../include;
cd -;

cd $sspfd_folder;
git pull;
make libsspfd.a;
cp libsspfd.a ../lib;
cp sspfd.h ../include;
cd -;

