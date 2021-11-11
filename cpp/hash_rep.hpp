/*
 Author: Shohei Kojima @ RIKEN
 This reads a relatively small fasta file and hash the sequences for both forward and reverse strands.
 */

#ifndef HASH_REP_HPP
#define HASH_REP_HPP

extern "C" {
namespace hasp_rep_hpp {

/*
 Reads a fasta file with multiple seqs and calc. hash values.
 */
int hash_rep(char* fa, int seed_len);

}  // namespace
}  // extern C
#endif
