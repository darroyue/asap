#ifndef RRR_BIT_VECTORS
#define RRR_BIT_VECTORS

#include <sdsl/bit_vectors.hpp>
#include "abstract_bv.h"

using namespace sdsl;

class RRRBV: public Bitvector {
  private:
    rrr_vector<> B;
    size_t vector_size;
    rrr_vector<>::rank_1_type *b_rank;
    rrr_vector<>::select_1_type *b_select;
  public:
    RRRBV (bit_vector input);
    RRRBV () {};
    ~RRRBV (void);
    int operator[](size_t i);
    int access (size_t i);
    int rank   (size_t i);
    int select (size_t i);
    int size ();
};

#endif
