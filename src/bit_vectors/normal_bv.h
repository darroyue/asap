#ifndef NORMAL_BIT_VECTORS
#define NORMAL_BIT_VECTORS

#include <sdsl/bit_vectors.hpp>
#include "abstract_bv.h"

using namespace sdsl;

class NBV: public Bitvector {
  private:
    bit_vector B;
    size_t vector_size;
    rank_support_v<1> *b_rank;
    bit_vector::select_1_type *b_select;
  public:
    NBV (bit_vector input);
    ~NBV (void);
    int operator[] (size_t i);
    int access (size_t i);
    int rank   (size_t i);
    int select (size_t i);
    int size ();
};

#endif
