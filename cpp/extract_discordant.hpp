/*
 Author: Shohei Kojima @ RIKEN
 This is just a header file of extract_discordant.cpp
 */

extern "C" {
namespace extract_discordant_hpp {

bool is_rep(std::string seq, uint64_t clip_len, const std::vector<uint32_t>& crepkmer, const ull& num_kmer);
int extract_discordant(std::string bam, std::string f_mainchr, std::string mk, std::string mi, std::string ref_fa,
                       std::string ourdir, int n_thread, bool is_cram);
int main(int argc, char* argv[]);

}
}
