#ifndef SD_BIT_VECTORS
#define SD_BIT_VECTORS

#include <sdsl/bit_vectors.hpp>
#include "abstract_bv.h"

using namespace sdsl;

class SDBV: public Bitvector {
  private:
    sd_vector<> B;
    size_t vector_size;
    sd_vector<>::rank_1_type *b_rank;
    sd_vector<>::select_1_type *b_select;
  public:
    SDBV (bit_vector input);
    ~SDBV (void);
    int operator[](size_t i);
    int access (size_t i);
    int rank   (size_t i);
    int select (size_t i);
    int size ();
};

#endif