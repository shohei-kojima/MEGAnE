/*
 Author: Shohei Kojima @ RIKEN
 This is just a header file of remove_multimapping_reads_from_fa.cpp
 */

extern "C" {

int remove_multimapping(char* in_mk, char* in_mi, char* in_fa, char* out_fa);
int main(int argc, char* argv[]);

}
