#include "iostream"

#include "TROOT.h"
#include "TVectorD.h"
#include "eventVGeneral.hpp"

using namespace std;


eventVGeneral::eventVGeneral(Double_t xmin, Double_t xmax,
                             Double_t ymin, Double_t ymax,
                             Double_t zmin, Double_t zmax): eventVolume(){
  // Double_t bottom_[] = {xmin, ymin, zmin};
  // bottom = TVectorD(3, bottom_);
  bottom[0] = xmin;
  bottom[1] = ymin;
  bottom[2] = zmin;
  // Double_t top_[] = {xmax, ymax, zmax};
  // top = TVectorD(3, top_);
  top[0] = xmax;
  top[1] = ymax;
  top[2] = zmax;
}

eventVGeneral::eventVGeneral(Double_t xlen,
                             Double_t ylen,
                             Double_t zlen): eventVolume() {
  // std::cout<< "Generate " <<xlen<< ", " << ylen << ", " << zlen << endl;
  eventVGeneral(-1*xlen/2, xlen/2,
                -1*ylen/2, ylen/2,
                -1*zlen/2, zlen/2);
}

Bool_t eventVGeneral::Inside(TVectorD position){
  Double_t x, y, z;
  x = position[0];
  y = position[1];
  z = position[2];

  if(bottom[0] > x || top[0] < x){
    return false;
  }
  if(bottom[1] > y || top[1] < y){
    return false;
  }
  if(bottom[2] > z || top[2] < z){
    return false;
  }
  return true;
}
