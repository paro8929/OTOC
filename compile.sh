g++ -fpermissive -fopenmp -std=c++11 -I. gen_matrices.cpp -o gen_mat
g++ -fpermissive -fopenmp -std=c++11 -I. -I/home/paulrom/numerics/eigen/eigen-3.3.7 evs.cpp -o evs
g++ -fpermissive -fopenmp -std=c++11 -I. otoc.cpp -o otoc
