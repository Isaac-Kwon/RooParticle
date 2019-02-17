#ifndef H_EVENTVOLUME
#define H_EVENTVOLUME 1;

#include "TROOT.h"
#include "TVectorT.h"
#include "TVectorD.h"

class eventVolume{
public:
  eventVolume();
  virtual Bool_t Inside(TVectorD position){return true;}
protected:
private:
};

#endif
