#ifndef H_RECORDER
#define H_RECORDER 1;

#include "TROOT.h"
#include "TString.h"
#include "iostream"

template<typename Element>
class recorder{
    public:
    recorder(){return;}
    recorder(Element * addr_){addr=addr_;}
    ~recorder(){}
    void    SetAddress(Element * addr_){ addr = addr_; }
    Element GetData(){return data;}
    void    Save(){ if(!addr){return;} Save(*addr); }
    void    Save(const Element data_){data=data_;}
    Bool_t  MaxSave();
    Bool_t  MaxSave(const Element data_);
    Bool_t  MinSave();
    Bool_t  MinSave(const Element data_);
    void    Print();
    protected:
    Bool_t   CMPSave(const Bool_t sgn);
    Bool_t   CMPSave(const Bool_t sgn, const Element data_);
    private:
    Element data;
    Element * addr = nullptr;
};

typedef recorder<Bool_t>   recorderB;
typedef recorder<Int_t>    recorderI;
typedef recorder<Long_t>   recorderL;
typedef recorder<Float_t>  recorderF;
typedef recorder<Double_t> recorderD;
#endif