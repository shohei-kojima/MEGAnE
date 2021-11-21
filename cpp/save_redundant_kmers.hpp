/*
 Author: Shohei Kojima @ RIKEN
 This is just a header file of save_redundant_kmers.cpp
 */

extern "C" {

int find_and_save_red_kmers(char* fa, char* f_fai, char* out_mk, char* out_mi);
int main(int argc, char* argv[]);

}
