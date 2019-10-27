#include "TROOT.h"
#include "event.hpp"

class oddevent : public event{
    public:
    void DeriveDT(Double_t dt = 1);
};