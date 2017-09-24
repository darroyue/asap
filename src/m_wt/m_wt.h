#ifndef M_WT_CLASS
#define M_WT_CLASS

#include <vector>
#include <unordered_map>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

template<class IntType = uint8_t>
class M_WT {
  private:
    sdsl::wt_ap<> m_wavelet_tree;
    std::unordered_map<IntType, unsigned> *m_map;
    std::vector<IntType> m_conversor;
    std::unordered_map<IntType, unsigned> m_rconversor;
    std::unordered_map<unsigned, bool> singleton;

  public:
    std::unordered_map<unsigned, unsigned> m_sizes;

    M_WT () {};
    M_WT (std::unordered_map<IntType, unsigned>, unsigned, unsigned);
    void m_default (std::unordered_map<IntType, unsigned>, unsigned);
    void m_by_freq (std::unordered_map<IntType, unsigned>, unsigned);
    void m_by_freq_simon (std::unordered_map<IntType, unsigned>, unsigned);
    void m_by_freq_rr (std::unordered_map<IntType, unsigned>, unsigned);

    unsigned map (IntType);
    IntType get_char_by_pos (unsigned);
    unsigned get_pos_by_char (IntType);

    unsigned access (unsigned);
    unsigned rank (unsigned, unsigned);
    unsigned select (unsigned, unsigned);
    unsigned long size ();

    bool is_singleton (unsigned);
    unsigned partitions ();
};

#include "m_wt.cpp"
#endif
