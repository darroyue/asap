#ifndef ABSTRACT_BIT_VECTOR
#define ABSTRACT_BIT_VECTOR
class Bitvector {
  public:
    virtual int operator[](size_t i) = 0;
    virtual int rank    (size_t i) = 0;
    virtual int select  (size_t i) = 0;
    virtual int access  (size_t i) = 0;
    virtual int size    () = 0;
    virtual ~Bitvector ( void ) = 0;
};
inline Bitvector::~Bitvector() { }
#endif