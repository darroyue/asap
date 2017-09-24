#include "arroyuelo_int.h"
#include "bit_vectors.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/vectors.hpp>

using namespace std;
using namespace sdsl;

template<class BitVectorClass, class IntType>
Arroyuelo<BitVectorClass, IntType>::Arroyuelo ( string input_file, unsigned method ) {
  IntType* s;
  // unordered_map<char, unsigned> freq;                        // freq[caracter] = frecuencia
  unordered_map<unsigned, unsigned> s_indexes;

  // Obtencion de las frequencias y de la cantidad de caracteres
  tie(s, text_length, freq) = Arroyuelo::readfile(input_file);
  m = new M_WT<IntType>(freq, text_length, method);

  unordered_map<unsigned, int_vector<>> s_vectors;
  unordered_map<unsigned, bit_vector> temp_bit_vectors;
  for ( auto &x: m->m_sizes ) {
    temp_bit_vectors[x.first] = bit_vector(text_length, 0);
    s_vectors[x.first] = int_vector<>(x.second);
  }

  // creacion de los bitvectors
  unsigned l;
  for ( unsigned i = 0; i < text_length; i++ ) {
    l = m->map(s[i]);
    temp_bit_vectors[l][i] = 1;
    s_vectors[l][s_indexes[l]++] = m->rank(m->get_pos_by_char(s[i]), l);
  }

  for ( auto &x: s_vectors )
    if ( ! m->is_singleton(x.first) )
      construct_im(s_wt_trees[x.first], x.second, 0);

  for ( auto &x: temp_bit_vectors )
    bit_vectors.insert(std::make_pair<unsigned, BitVectorClass*>((unsigned int)x.first, new BitVectorClass(x.second)));
}

template<class BitVectorClass, class IntType>
Arroyuelo<BitVectorClass, IntType>::~Arroyuelo ( void ) {
}

template<class BitVectorClass, class IntType>
IntType Arroyuelo<BitVectorClass, IntType>::access ( unsigned position ) {
  unsigned l = 0;
  for ( auto &x: bit_vectors )
    if ( (x.second)->access(position) == 1 ) {
      l = x.first;
      break;
    }

  if ( m->is_singleton(l) ) return m->get_char_by_pos(m->select(1, l));

  unsigned k = (bit_vectors[l])->rank(position);
  return m->get_char_by_pos(m->select(s_wt_trees[l][k] + 1, l));
}

template<class BitVectorClass, class IntType>
unsigned Arroyuelo<BitVectorClass, IntType>::rank ( IntType target, unsigned position ) {
  int l = m->map(target);
  // unsigned k = (bit_vectors[l])->rank(position + 1);
  unsigned k = (bit_vectors[l])->rank(position);

  if ( m->is_singleton(l) ) return k;

  unsigned c = m->rank(m->get_pos_by_char(target), l);
  return s_wt_trees[l].rank(k, c);
}

template<class BitVectorClass, class IntType>
int Arroyuelo<BitVectorClass, IntType>::select ( IntType target, unsigned index) {
  int l = m->map(target);
  if ( l == -1 ) return -1;

  if ( m->is_singleton(l) ) return (bit_vectors[l])->select(index);

  unsigned c = m->rank(m->get_pos_by_char(target), l);
  unsigned k = s_wt_trees[l].select(index, c) + 1;
  return (bit_vectors[l])->select(k);
}

template<class BitVectorClass, class IntType>
IntType* Arroyuelo<BitVectorClass, IntType>::waccess ( unsigned start, unsigned end ) {
  if ( end < start ) {
    unsigned temp = start;
    start = end;
    end = temp;
  }
  end = end < text_length ? end : text_length - 1;
  IntType* response = new IntType[end - start + 1];
  unsigned length = end - start + 1;
  unsigned relevant_bits, before_start_rank, end_rank;
  for ( auto &x: bit_vectors ) {
    before_start_rank = start == 0 ? 0 : (x.second)->rank(start);
    end_rank = (x.second)->rank(end+1);
    relevant_bits = end_rank - before_start_rank;

    if ( relevant_bits > 0 ) {
      unsigned nextOne, l = x.first, k = before_start_rank;
      for ( unsigned i = 1; nextOne = (x.second)->select(before_start_rank + i), nextOne <= end && relevant_bits > 0 ; i++, relevant_bits-- ) {
        response[nextOne - start] = m->get_char_by_pos(m->select( (m->is_singleton(l) ? 0 : s_wt_trees[l][k + 1]), l));
        //response[nextOne - start] = m->get_char_by_pos(m->select( (m->is_singleton(l) ? 0 : s_wt_trees[l][k]) + 1, l));
        length--;
        k += 1;
      }
    }
    if ( length == 0 ) break;
  }

  return response;
}

template<class BitVectorClass, class IntType>
unsigned long Arroyuelo<BitVectorClass, IntType>::size () {
  unsigned long sum = m->size();
  for ( auto &x: s_wt_trees )
    if ( x.second.sigma > 1 )
      sum += size_in_bytes(x.second);
  for ( auto &x: bit_vectors )
    sum += (x.second)->size();
  return sum;
}

template<class BitVectorClass, class IntType>
unsigned long Arroyuelo<BitVectorClass, IntType>::m_size () {
  return m->size();
}

template<class BitVectorClass, class IntType>
unsigned Arroyuelo<BitVectorClass, IntType>::partitions () {
  return m->partitions();
}

template<class BitVectorClass, class IntType>
unsigned long Arroyuelo<BitVectorClass, IntType>::b_size () {
  unsigned long sum = 0;
  for ( auto &x: bit_vectors )
    sum += (x.second)->size();
  return sum;
}

template<class BitVectorClass, class IntType>
unsigned long Arroyuelo<BitVectorClass, IntType>::sl_size () {
  unsigned long sum = 0;
  for ( auto &x: s_wt_trees )
    if ( x.second.sigma > 1 )
      sum += size_in_bytes(x.second);
  return sum;
}

template<class BitVectorClass, class IntType>
tuple<IntType*, unsigned, unordered_map<IntType, unsigned>> Arroyuelo<BitVectorClass, IntType>::readfile ( string input_file ) {
  FILE *file;
  IntType* buffer;
  unsigned length, result;
  unordered_map<IntType, unsigned> freq;

  file = fopen(input_file.c_str(), "rb");
  fseek(file, 0, SEEK_END);
  length = ftell(file)/sizeof(IntType);
  rewind(file);
  buffer = new IntType[length];
  result = fread(buffer, sizeof(IntType), length, file);

  if ( result != length ) { cout << "fread error" << endl; exit(3);}

  fclose(file);

  for ( unsigned i = 0; i < length; i++ )
    freq[buffer[i]]++;

  return make_tuple (buffer, length, freq);
}

template<class BitVectorClass, class IntType>
unsigned Arroyuelo<BitVectorClass, IntType>::length () {
  return text_length;
}

template<class BitVectorClass, class IntType>
unsigned Arroyuelo<BitVectorClass, IntType>::alphabet_size() {
  return freq.size();
}

template<class BitVectorClass, class IntType>
std::unordered_map<IntType, unsigned> Arroyuelo<BitVectorClass, IntType>::frequency() {
  return freq;
}

template<class BitVectorClass, class IntType>
vector<IntType> Arroyuelo<BitVectorClass, IntType>::alphabet () {
  vector<IntType> alphabet;
  for ( auto &x: freq )
    alphabet.push_back(x.first);
  return alphabet;
}

template<class BitVectorClass, class IntType>
std::tuple<IntType, double, bool, double, double, double, double> Arroyuelo<BitVectorClass, IntType>::access_timecheck ( unsigned position ) {
  double l_find_t, singleton_t, bv_rank_t, wt_access_t, m_t;
  clock_t begin, end;

  begin = clock();
  unsigned l;
  for ( auto &x: bit_vectors )
    if ( (x.second)->access(position) == 1 ) {
      l = x.first;
      break;
    }
  end = clock();
  l_find_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  if ( m->is_singleton(l) ) {
    begin = clock();
    IntType answer = m->get_char_by_pos(m->select(1, l));
    end = clock();
    m_t = double(end-begin) / CLOCKS_PER_SEC;

    return make_tuple(answer, l_find_t, true, 0, 0, 0, m_t);
  }
  end = clock();
  singleton_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned k = (bit_vectors[l])->rank(position);
  end = clock();
  bv_rank_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  auto temp = s_wt_trees[l][k];
  end = clock();
  wt_access_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  IntType temp2 = m->get_char_by_pos(m->select(temp + 1, l));
  end = clock();
  m_t = double(end-begin) / CLOCKS_PER_SEC;
  return make_tuple(temp2, l_find_t, false, singleton_t, bv_rank_t, wt_access_t, m_t);
}

template<class BitVectorClass, class IntType>
std::tuple<unsigned, double, double, bool, double, double, double> Arroyuelo<BitVectorClass, IntType>::rank_timecheck ( IntType target, unsigned position ) {
  double l_find_t, bv_rank_t, singleton_t, m_rank_t, wt_rank_t;
  clock_t begin, end;

  begin = clock();
  int l = m->map(target);
  end = clock();
  l_find_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned k = (bit_vectors[l])->rank(position + 1);
  end = clock();
  bv_rank_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  if ( m->is_singleton(l) ) {
    return make_tuple(k, l_find_t, bv_rank_t, true, 0, 0, 0);
  }
  end = clock();
  singleton_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned c = m->rank(m->get_pos_by_char(target), l);
  end = clock();
  m_rank_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned temp = s_wt_trees[l].rank(k, c);
  end = clock();
  wt_rank_t = double(end-begin) / CLOCKS_PER_SEC;

  return make_tuple(temp, l_find_t, bv_rank_t, false, singleton_t, m_rank_t, wt_rank_t);
}

template<class BitVectorClass, class IntType>
std::tuple<int, double, bool, double, double, double, double> Arroyuelo<BitVectorClass, IntType>::select_timecheck ( IntType target, unsigned index) {
  double find_l_t, singleton_t, m_rank_t, wt_select_t, bv_select_t;
  clock_t begin, end;

  begin = clock();
  int l = m->map(target);
  end = clock();
  find_l_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  if ( m->is_singleton(l) ) {
    begin = clock();
    int answer = (bit_vectors[l])->select(index);
    end = clock();
    bv_select_t = double(end-begin) / CLOCKS_PER_SEC;
    return make_tuple(answer, find_l_t, true, 0, 0, 0, bv_select_t);
  }
  end = clock();
  singleton_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned c = m->rank(m->get_pos_by_char(target), l);
  end = clock();
  m_rank_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned k = s_wt_trees[l].select(index, c) + 1;
  end = clock();
  wt_select_t = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  int temp = (bit_vectors[l])->select(k);
  end = clock();
  bv_select_t = double(end-begin) / CLOCKS_PER_SEC;

  return make_tuple(temp, find_l_t, false, singleton_t, m_rank_t, wt_select_t, bv_select_t);
}


template<class BitVectorClass, class IntType>
std::tuple<IntType, bool, unsigned, double, double> Arroyuelo<BitVectorClass, IntType>::accessTime ( unsigned position ) {
  IntType answer;
  clock_t begin = clock(), end;
  unsigned l = 0;
  for ( auto &x: bit_vectors )
    if ( (x.second)->access(position) == 1 ) {
      l = x.first;
      break;
    }
  end = clock();

  double b_and_s_time = double(end-begin) / CLOCKS_PER_SEC, m_time;

  begin = clock();
  if ( m->is_singleton(l) ) {
    answer = m->get_char_by_pos(m->select(1, l));
    end = clock();

    m_time = double(end-begin) / CLOCKS_PER_SEC;
    return std::make_tuple(answer, true, l, m_time, b_and_s_time);
  }
  end = clock();

  m_time = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned k = (bit_vectors[l])->rank(position);
  end = clock();

  b_and_s_time += double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  answer = m->get_char_by_pos(m->select(s_wt_trees[l][k] + 1, l));
  end = clock();

  m_time += double(end-begin) / CLOCKS_PER_SEC;
  return std::make_tuple(answer, true, l, m_time, b_and_s_time);
}

template<class BitVectorClass, class IntType>
std::tuple<unsigned, bool, unsigned, double, double> Arroyuelo<BitVectorClass, IntType>::rankTime ( IntType target, unsigned position ) {
  unsigned answer;
  clock_t begin = clock(), end;
  int l = m->map(target);
  end = clock();

  double m_time = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned k = (bit_vectors[l])->rank(position + 1);
  end = clock();

  double b_and_s_time = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  if ( m->is_singleton(l) ) {
    answer = k;
    end = clock();

    m_time += double(end-begin) / CLOCKS_PER_SEC;
    return std::make_tuple(answer, true, l, m_time, b_and_s_time);
  }

  unsigned c = m->rank(m->get_pos_by_char(target), l);
  end = clock();

  m_time += double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  answer = s_wt_trees[l].rank(k, c);
  end = clock();

  b_and_s_time += double(end-begin) / CLOCKS_PER_SEC;
  return std::make_tuple(answer, true, l, m_time, b_and_s_time);
}

template<class BitVectorClass, class IntType>
std::tuple<int, bool, unsigned, double, double> Arroyuelo<BitVectorClass, IntType>::selectTime ( IntType target, unsigned index) {
  unsigned answer;
  double m_time, b_and_s_time;
  clock_t begin = clock(), end;
  int l = m->map(target);
  if ( l == -1 ) {
    return std::make_tuple(-1, false, -1, 0, 0);
  }

  if ( m->is_singleton(l) ) {
    end = clock();
    m_time = double(end-begin) / CLOCKS_PER_SEC;

    begin = clock();
    answer = (bit_vectors[l])->select(index);
    end = clock();

    b_and_s_time = double(end-begin) / CLOCKS_PER_SEC;
    return std::make_tuple(answer, true, l, m_time, b_and_s_time);
  }

  unsigned c = m->rank(m->get_pos_by_char(target), l);
  end = clock();

  m_time = double(end-begin) / CLOCKS_PER_SEC;

  begin = clock();
  unsigned k = s_wt_trees[l].select(index, c) + 1;
  answer = (bit_vectors[l])->select(k);
  end = clock();

  b_and_s_time = double(end-begin) / CLOCKS_PER_SEC;
  return std::make_tuple(answer, true, l, m_time, b_and_s_time);
}

template<class BitVectorClass, class IntType>
std::tuple<IntType*, unordered_map<unsigned, double>> Arroyuelo<BitVectorClass, IntType>::waccessTime ( unsigned start, unsigned end ) {
  double initialTime;
  unordered_map<unsigned, double> partitionTime;
  clock_t begin = clock(), clockEnd;

  if ( end < start ) {
    unsigned temp = start;
    start = end;
    end = temp;
  }
  end = end < text_length ? end : text_length - 1;
  IntType* response = new IntType[end - start + 1];
  unsigned length = end - start + 1;
  unsigned relevant_bits, before_start_rank, end_rank;

  clockEnd = clock();

  initialTime = double(clockEnd-begin) / CLOCKS_PER_SEC;

  for ( auto &x: bit_vectors ) {
    begin = clock();

    before_start_rank = start == 0 ? 0 : (x.second)->rank(start);
    end_rank = (x.second)->rank(end+1);
    relevant_bits = end_rank - before_start_rank;

    if ( relevant_bits > 0 ) {
      unsigned nextOne, l = x.first, k = before_start_rank;
      for ( unsigned i = 1; nextOne = (x.second)->select(before_start_rank + i), nextOne <= end && relevant_bits > 0 ; i++, relevant_bits-- ) {
        response[nextOne - start] = m->get_char_by_pos(m->select( (m->is_singleton(l) ? 0 : s_wt_trees[l][k]) + 1, l));
        length--;
        k += 1;
      }
    }

    clockEnd = clock();

    partitionTime[x.first] = double(clockEnd-begin) / CLOCKS_PER_SEC;
    if ( length == 0 ) break;
  }

  for ( auto &x: partitionTime )
    x.second += initialTime;

  return std::make_tuple(response, partitionTime);
}
