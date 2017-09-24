#include "normal_bv.h"

NBV::NBV (bit_vector input) {
  vector_size = input.size();
  B = input;
  b_rank = new rank_support_v<1>(&B);
  b_select = new bit_vector::select_1_type(&B);
}

NBV::~NBV (void) {
  delete b_rank;
  delete b_select;
}

int NBV::operator[](size_t i) {
  return i < vector_size ? B[i] : 0;
}

int NBV::access (size_t i) {
  return i < vector_size ? B[i] : 0;
}

int NBV::rank   (size_t i) {
  return (*b_rank)(i);
}

int NBV::select (size_t i) {
  return i < vector_size ? (*b_select)(i) : 0;
}

int NBV::size () {
  return size_in_bytes(*b_rank) + size_in_bytes(*b_select) + size_in_bytes(B);
}