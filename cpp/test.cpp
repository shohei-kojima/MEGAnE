#include <iostream>
#include "hash_rep.hpp"
using namespace hasp_rep_hpp;


int main(int argc, char* argv[]) {
    if (argc <= 2) {
        std::cerr << "Please specify a fasta file and seq length to be hashed." << std::endl;
        return 1;
    }
    char* fa=argv[1];
    int seed_len=atoi(argv[2]);
    std::cout << "Reading " << argv[1] << " ..." << std::endl;
    
    int ret=hash_rep(fa, seed_len);
    return 0;
}


/*
 g++ -c hash_rep.cpp
 g++ -c test.cpp
 g++ -o test hash_rep.o test.o
 g++ -o test test.cpp hash_rep.cpp
 -L/home/kooojiii/results/2021/prog_develop/MEGAnE/cpp -I/home/kooojiii/results/2021/prog_develop/MEGAnE/cpp
 /home/kooojiii/Documents/prog_develop/singularity/MEGAnE/210118_1/snd/usr/local/bin/MEGAnE/lib/humrepsub.fa 11
 */
