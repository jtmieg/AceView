nvcc -x c++ -g -O4  -DACEDB4 \
 -I/usr/include/c++/8 -I/usr/include/c++/8/x86_64-redhat-linux \
 -I. -I.. -I../wh -I../wsa -I/netopt/sge/include \
 -DAC_TEST -DOPTERON -D_FILE_OFFSET_BITS=64 -c -Wno-deprecated-gpu-targets  -o sa.main.o sa.main.c
