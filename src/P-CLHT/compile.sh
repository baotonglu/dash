BASE=$1

cd external/ssmem
make clean
make
cp libssmem.a ../lib/

cd ../../
make clean
make libclht_${BASE}_res.a
