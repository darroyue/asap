#include "sd_bv.h"

SDBV::SDBV (bit_vector input) {
  vector_size = input.size();
  B = sd_vector<>(input);
  b_rank = new sd_vector<>::rank_1_type(&B);
  b_select = new sd_vector<>::select_1_type(&B);
}

SDBV::~SDBV (void) {
  delete b_rank;
  delete b_select;
}

int SDBV::operator[](size_t i) {
  return i < vector_size ? B[i] : 0;
}

int SDBV::access (size_t i) {
  return i < vector_size ? B[i] : 0;
}

int SDBV::rank (size_t i) {
  return (*b_rank)(i);
}

int SDBV::select (size_t i) {
  return i < vector_size ? (*b_select)(i) : 0;
}

int SDBV::size () {
  return size_in_bytes(*b_rank) + size_in_bytes(*b_select) + size_in_bytes(B);
}