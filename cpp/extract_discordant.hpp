/*
 Author: Shohei Kojima @ RIKEN
 This is just a header file of extract_discordant.cpp
 */

extern "C" {

int extract_discordant(std::string bam, std::string f_mainchr, std::string mk, std::string mi, std::string ref_fa,
                       std::string ourdir, int n_thread, bool is_cram);
int main(int argc, char* argv[]);

}
