make -f Makefile.cvs
mkdir release
cd release
CXXFLAGS="-O3 -funroll-loops -mfpmath=sse -msse -msse2" LDFLAGS="-lpthread" ../configure
make -j1

