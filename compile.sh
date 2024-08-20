g++ -std=c++11 MolSuperposition.cpp -I/usr/local/include/openbabel3 \
    -L/usr/local/lib \
    -I/usr/include/eigen3 \
    -Wl,-rpath=/usr/local/lib -l openbabel -o molsuperpos


