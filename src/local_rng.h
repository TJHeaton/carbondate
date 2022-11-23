#ifndef LOCALRNG_H
#define LOCALRNG_H

#include "R.h"

// Get and set the random number generator set on entering and exiting functions
// Important so we get the same stream of random numbers when interchanging
// between R and C++ as we would by running the code only in R.
class local_rng {
public:
  local_rng() {
    GetRNGstate();
  }

  ~local_rng(){
    PutRNGstate();
  }
};

#endif
