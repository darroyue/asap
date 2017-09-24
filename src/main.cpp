#include "arroyuelo_int.h"
#include "bit_vectors.h"

int main(int argc, char *argv[]) {
  // Constructor:
  // templateArg0: Bitvector to use
  // templateArg1: type of encoding used on source file's words
  // arg0: path to encoded words file
  // arg1: partition method; 0: sparse, 1: dense, 2: dense-l_min
  Arroyuelo<SDBV, uint32_t> AWT("path_to_file", 0);
  return 0;
}
