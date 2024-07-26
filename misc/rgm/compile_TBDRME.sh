CC -std=c++14 -I/global/homes/j/jherko/Research/wigxjpf-1.11/inc -I/global/common/software/m2032/shared/spack/install/cray-cnl7-haswell/boost-1.78.0-fats5tl/include -I/global/homes/j/jherko/Research/lib/eigen-3.4.0 -O3  -DNDEBUG  -DDTU3R3_CACHE  -DHAVE_80BIT_LONG_DOUBLE -DHAVE_FLOAT128 -I/global/homes/j/jherko/Research/lib/SU3lib/include -c -o TBDRME.o TBDRME.cpp

CC -o TBDRME TBDRME.o -L/global/homes/j/jherko/Research/lib/SU3lib/src -lSU3 -L/global/homes/j/jherko/Research/wigxjpf-1.11/lib -lwigxjpf -I/global/homes/j/jherko/Research/lib/eigen-3.4.0
