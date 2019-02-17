#ifndef H_EVENTVGENERAL
#define H_EVENTVGENERAL 1;

#include "TROOT.h"
#include "TVectorD.h"
#include "eventVolume.hpp"

class eventVGeneral: public eventVolume{
public:
  eventVGeneral(Double_t xmin, Double_t xmax,
                Double_t ymin, Double_t ymax,
                Double_t zmin, Double_t zmax);
  eventVGeneral(Double_t xlen,
                Double_t ylen,
                Double_t zlen);
  Bool_t Inside(TVectorD position);
protected:
private:
  TVectorD bottom = TVectorD(3);
  TVectorD top = TVectorD(3);
};

#endif
