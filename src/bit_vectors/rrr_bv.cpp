#include "rrr_bv.h"

RRRBV::RRRBV (bit_vector input) {
  vector_size = input.size();
  B = rrr_vector<>(input);
  b_rank = new rrr_vector<>::rank_1_type(&B);
  b_select = new rrr_vector<>::select_1_type(&B);
}

RRRBV::~RRRBV (void) {
  // delete b_rank;
  // delete b_select;
}

int RRRBV::operator[](size_t i) {
  return i < vector_size ? B[i] : 0;
}

int RRRBV::access (size_t i) {
  return i < vector_size ? B[i] : 0;
}

int RRRBV::rank   (size_t i) {
  return (*b_rank)(i);
}

int RRRBV::select (size_t i) {
  return i < vector_size ? (*b_select)(i) : 0;
}

int RRRBV::size () {
  return size_in_bytes(*b_rank) + size_in_bytes(*b_select) + size_in_bytes(B);
}