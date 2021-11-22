/*
 Author: Shohei Kojima @ RIKEN
 This is just a header file of remove_multimapping_reads_from_fa.cpp
 */

extern "C" {

int find_and_save_all_kmers(char* fa, char* out_mk, char* out_mi);
int main(int argc, char* argv[]);

}
