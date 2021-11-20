
/*
 Author: Shohei Kojima @ RIKEN
 This is a minimal fasta index parser.
 */

#ifndef PARSE_FAI_HPP
#define PARSE_FAI_HPP

#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cassert>

typedef unsigned long long ull;

extern "C" {
namespace parse_fai_hpp {

/*
 this holds one line of fai file
 */
struct fai_info {
    int chr_name_len     = 0;
    char* chr_name       = nullptr;
    ull chr_len          = 0;
    ull seq_start        = 0;
    ull seq_len_per_line = 0;
    ull byte_per_line    = 0;
    
    fai_info(std::vector<std::string> v) {
        assert(v.size() == 5);
        int hlen               = v[0].length();
        this->chr_name_len     = hlen;
        this->chr_name         = new char[hlen + 1];
        std::strcpy(this->chr_name, v[0].c_str());
        this->chr_len          = std::stoull(v[1]);
        this->seq_start        = std::stoull(v[2]);
        this->seq_len_per_line = std::stoull(v[3]);
        this->byte_per_line    = std::stoull(v[4]);
    };
};

/*
 usage:
    char* f_fai="file.fa.fai";
    fai_parser fai;
    int chr_num=fai.load_fai(f_fai);
    fai_info* chr_info=fai.get_chr_info(0);  // chr1
    cout << chr_info->chr_name_len << endl;
    cout << chr_info->chr_name << endl;
    cout << chr_info->chr_len << endl;
    cout << chr_info->seq_start << endl;
    cout << chr_info->seq_len_per_line << endl;
    cout << chr_info->byte_per_line << endl;
 */
class fai_parser {
    std::vector<fai_info*> chr_infos;
public:
    int load_fai(char* f);
    int get_chr_num();
    fai_info* get_chr_info(int n);
};

int fai_parser::load_fai(char* f) {
    // open fai
    std::ifstream infile;
    infile.open(f);
    if (! infile.is_open()) {
        return 1;
    }
    
    // load fai
    std::string line;
    int n=0;
    while (std::getline(infile, line)) {
        std::stringstream ss{line};
        std::string buf;
        std::vector<std::string> v;
        while (std::getline(ss, buf, '\t')) {
          v.push_back(buf);
        }
        fai_info* info= new fai_info(v);
        this->chr_infos.push_back(info);
        n++;
    }
    
    // close
    infile.close();
    
    return n;
}

int fai_parser::get_chr_num() {
    return this->chr_infos.size();
}

fai_info* fai_parser::get_chr_info(int n) {
    return this->chr_infos[n];
}

}  // namespace
}  // extern C
#endif
