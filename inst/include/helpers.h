#ifndef __HELPERS__
#define __HELPERS__

//===============================================//
//
// Contains a variety of global helper functions
//
//===============================================//

inline double max(double x, double y) {
  // Args:
  //  * x, y: Scalars.
  // Returns:
  //  * max(x, y).
  if (x >= y) {
    return x;
  } else {
    return y;
  }
}

#endif
