#include "recorder.hpp"




template <typename Element>
void recorder<Element>::Print(){
    TString prt;
    if(addr == nullptr){
        std::cout<<"recorderT::Print : Target variable address is not saved."<<std::endl;
        prt = "nullptr";
    }else{
        prt = *addr;
    }

    std::cout << "RecorderT | data_addr: " << prt << " | data_saved " << data << std::endl;
    return;
}

template <typename Element> 
Bool_t recorder<Element>::CMPSave(const Bool_t sgn, const Element data_){
// if sgn kTRUE, save higher value. if kFALSE, save lower value.
    if((sgn && data_>data) || (sgn && data_<data)){
        Save(data_);
        return kTRUE; //If saved, return kTRUE
    }
    return kFALSE; //If NOT saved, return kFALSE
}

template <typename Element> 
Bool_t recorder<Element>::CMPSave(const Bool_t sgn){
    return CMPSave(sgn, *addr);
}

template <typename Element> 
Bool_t recorder<Element>::MaxSave(const Element data_){
    return CMPSave(kTRUE, data_);
}

template <typename Element> 
Bool_t recorder<Element>::MaxSave(){
    return CMPSave(kTRUE);
}

template <typename Element> 
Bool_t recorder<Element>::MinSave(const Element data_){
    return CMPSave(kFALSE, data_);
}

template <typename Element> 
Bool_t recorder<Element>::MinSave(){
    return CMPSave(kFALSE);
}


template class  recorder<Bool_t>;  
template class  recorder<Int_t>;
template class  recorder<Long_t>;  
template class  recorder<Float_t>; 
template class  recorder<Double_t>; 

   


