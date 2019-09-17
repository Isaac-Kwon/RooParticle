#ifndef H_EVENTVOLUME
#define H_EVENTVOLUME 1;

#include "TROOT.h"
#include "TVector3.h"

class eventVolume{
public:
  eventVolume();
  virtual Bool_t Inside(TVector3 position){return true;}
protected:
private:
};

class eventVGeneral: public eventVolume{
public:
  eventVGeneral(Double_t xmin, Double_t xmax,
                Double_t ymin, Double_t ymax,
                Double_t zmin, Double_t zmax);
  eventVGeneral(Double_t xlen, Double_t ylen, Double_t zlen);
  Bool_t Inside(TVector3 position);
protected:
private:
  TVector3 bottom;
  TVector3 top;
};

#endif
