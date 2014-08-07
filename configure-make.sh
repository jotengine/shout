make -f Makefile.cvs
mkdir release
cd release
CXXFLAGS="-O3 -funroll-loops -mfpmath=sse -msse -msse2" LDFLAGS="-lpthread" CXX="/usr/local/bin/g++-4.8" ../configure
make -j1