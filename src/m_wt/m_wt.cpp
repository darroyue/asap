#include "m_wt.h"

using namespace std;
using namespace sdsl;

template<class IntType>
M_WT<IntType>::M_WT (unordered_map<IntType, unsigned> freq, unsigned n, unsigned method) {
  m_map = new std::unordered_map<IntType, unsigned>();
  switch ( method ) {
    // using buckets with ceiling of log... formula
    case 0:
      M_WT::m_default(freq, n);
      break;
    // by frequency rank
    case 1:
      M_WT::m_by_freq(freq, n);
      break;
    // by frequency rank with Simon Gog's adjustement
    case 2:
      M_WT::m_by_freq_simon(freq, n);
      break;
    // by frequency rank && round robin
    case 3:
      M_WT::m_by_freq_rr(freq, n);
      break;
  }
}

template<class IntType>
void M_WT<IntType>::m_default (unordered_map<IntType, unsigned> freq, unsigned n) {
  int_vector<> m_vector(freq.size());
  double logn = log2(n);
  unsigned i = 0;
  for ( auto &val: freq) {
    unsigned ceiling = ceil( (logn - log2(val.second)) * logn );
    if ( singleton.count(ceiling)  == 0 )
      singleton[ceiling] = true;
    else
      singleton[ceiling] = false;
    m_vector[i] = ceiling;
    (*m_map)[val.first] = ceiling;
    m_conversor.push_back(val.first);
    m_rconversor[val.first] = i++;
    m_sizes[ceiling] += val.second;
  }

  construct_im(m_wavelet_tree, m_vector);
}

template<class IntType>
void M_WT<IntType>::m_by_freq (unordered_map<IntType, unsigned> freq, unsigned n) {
  int_vector<> m_vector(freq.size());

  vector<pair<IntType, unsigned>> byFreq;
  for ( auto &val: freq )
    byFreq.push_back(make_pair(val.first, val.second));

  sort(byFreq.begin(), byFreq.end(), [](const pair<IntType, unsigned> &left, const pair<IntType, unsigned> &right) {
    return left.second > right.second;
  });

  unsigned currentSize = 1, count = 0, i = 0;
  for ( auto &val: byFreq ) {
    double l = log2(currentSize);
    if ( singleton.count(l)  == 0 )
      singleton[l] = true;
    else
      singleton[l] = false;
    m_vector[i] = l;
    (*m_map)[val.first] = l;
    m_conversor.push_back(val.first);
    m_rconversor[val.first] = i++;
    m_sizes[l] += val.second;

    if ( ++count == currentSize ) {
      count = 0;
      currentSize *= 2;
    }
  }

  construct_im(m_wavelet_tree, m_vector);
}

template<class IntType>
void M_WT<IntType>::m_by_freq_simon (unordered_map<IntType, unsigned> freq, unsigned n) {
  int_vector<> m_vector(freq.size());

  vector<pair<IntType, unsigned>> byFreq;
  for ( auto &val: freq )
    byFreq.push_back(make_pair(val.first, val.second));

  sort(byFreq.begin(), byFreq.end(), [](const pair<IntType, unsigned> &left, const pair<IntType, unsigned> &right) {
    return left.second > right.second;
  });

  unsigned freq_size = freq.size();
  unsigned m_singleton_class_cnt = log2(freq_size);

  for ( unsigned i = 0; i < m_singleton_class_cnt; i++ ) {
    singleton[i] = true;
    m_vector[i] = i;
    (*m_map)[byFreq[i].first] = i;
    m_conversor.push_back(byFreq[i].first);
    m_rconversor[byFreq[i].first] = i;
    m_sizes[i] = byFreq[i].second;
  }

  unsigned currentSize = 2, count = 0;
  unsigned currentL = m_singleton_class_cnt;
  for ( unsigned i = m_singleton_class_cnt; i < freq_size; i++ ) {
    singleton[currentL] = false;
    m_vector[i] = currentL;
    (*m_map)[byFreq[i].first] = currentL;
    m_conversor.push_back(byFreq[i].first);
    m_rconversor[byFreq[i].first] = i;
    m_sizes[currentL] += byFreq[i].second;

    if ( ++count == currentSize ) {
      count = 0;
      currentSize *= 2;
      currentL++;
    }
  }

  construct_im(m_wavelet_tree, m_vector);
}

template<class IntType>
void M_WT<IntType>::m_by_freq_rr (unordered_map<IntType, unsigned> freq, unsigned n) {
  int_vector<sizeof(IntType)*8> m_vector(freq.size());

  vector<pair<IntType, unsigned>> byFreq;
  for ( auto &val: freq )
    byFreq.push_back(make_pair(val.first, val.second));

  sort(byFreq.begin(), byFreq.end(), [](const pair<IntType, unsigned> &left, const pair<IntType, unsigned> &right) {
    return left.second < right.second;
  });

  /*WIP*/
}

template<class IntType>
unsigned M_WT<IntType>::map (IntType c) {
  return (*m_map)[c];
}

template<class IntType>
IntType M_WT<IntType>::get_char_by_pos (unsigned pos) {
  return m_conversor[pos];
}

template<class IntType>
unsigned M_WT<IntType>::get_pos_by_char (IntType c) {
  return m_rconversor[c];
}

template<class IntType>
unsigned M_WT<IntType>::access (unsigned pos) {
  return m_wavelet_tree[pos];
}

template<class IntType>
unsigned M_WT<IntType>::rank (unsigned pos, unsigned  l) {
  return m_wavelet_tree.rank(pos, l);
}

template<class IntType>
unsigned M_WT<IntType>::select (unsigned target, unsigned  l) {
  return m_wavelet_tree.select(target, l);
}

template<class IntType>
unsigned long M_WT<IntType>::size () {
  unsigned long sigma = m_map->size();
  unsigned long size_of_m_map = sigma * ( 4 + sizeof(IntType) );
  unsigned long size_of_m_conversor = sigma * sizeof(IntType);
  unsigned long size_of_m_rconversor = sigma * ( 4 + sizeof(IntType) );

  return size_in_bytes(m_wavelet_tree) + size_of_m_map + size_of_m_conversor + size_of_m_rconversor;
}

template<class IntType>
bool M_WT<IntType>::is_singleton(unsigned l) {
  return singleton[l];
}

template<class IntType>
unsigned M_WT<IntType>::partitions () {
  return singleton.size();
}