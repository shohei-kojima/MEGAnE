/*
 Author: Shohei Kojima @ RIKEN
 This is a minimum prog to generate complementary DNA seq.
 This is a header-only lib.
 */

#ifndef COMPLEMENTARY_SEQ_HPP
#define COMPLEMENTARY_SEQ_HPP

extern "C" {
namespace complementary_seq_hpp {

typedef unsigned long long ull;

static char complement[128]={'N'};

void comp_init() {
    complement['A']='T'; complement['a']='T';
    complement['C']='G'; complement['c']='G';
    complement['G']='C'; complement['g']='C';
    complement['T']='A'; complement['t']='A';
    complement['U']='A'; complement['u']='A';
    complement['M']='K'; complement['m']='K';
    complement['R']='Y'; complement['r']='Y';
    complement['W']='W'; complement['w']='W';
    complement['S']='S'; complement['s']='S';
    complement['Y']='R'; complement['y']='R';
    complement['K']='M'; complement['k']='M';
    complement['V']='B'; complement['v']='B';
    complement['H']='D'; complement['h']='D';
    complement['D']='H'; complement['d']='H';
    complement['B']='V'; complement['b']='V';
    complement['N']='N'; complement['n']='N';
}


/*
 Takes a sequence and returns complementary seq
 Args:
   1) pointer to sequence
   2) length of sequence
 */
char* complementary_seq(const char* seq, ull& leng) {
    char* cseq= new char[leng + 1];
    for (ull i=0; i < leng; i++) {
        cseq[i]=complement[seq[leng - i - 1]];
    }
    cseq[leng]='\0';
    return cseq;
}


}  // namespace
}  // extern C
#endif
