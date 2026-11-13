/// \file recorder.hpp
/// \brief Small value holder that remembers one number, optionally read from a variable (\ref recorder).

#ifndef H_RECORDER
#define H_RECORDER 1;

#include "TROOT.h"
#include "TString.h"
#include "iostream"

/// \brief Remembers one value of type \c Element, optionally pulled from a watched variable.
///
/// Intended use: keep the extreme of a quantity over a run (for example the
/// minimum distance of closest approach). The member functions are defined in
/// recorder.cpp and instantiated there only for the typedefs recorderB (Bool_t),
/// recorderI (Int_t), recorderL (Long_t), recorderF (Float_t) and recorderD
/// (Double_t); other element types fail to link. All declared members are
/// defined. The recorder does not own the watched variable.
///
/// For a minimum search start the recorder with the largest possible value
/// (std::numeric_limits<Double_t>::max()) and feed it with MinSave(), as
/// test/main_3body_oneline.cpp and test/main_rutherford_recorder.cpp do.
///
/// \note Until 2026-11 a typo in the comparison made MinSave() never save and
///       MaxSave() save any differing value; output written before that holds
///       the initial value instead of the minimum.
template<typename Element>
class recorder{
    public:
    /// \brief Empty recorder: no watched variable, stored value uninitialised (for built-in types).
    recorder(){return;}
    /// \brief Recorder with an initial stored value.
    /// \param data_ initial value (for a minimum search use the largest possible value)
    recorder(Element   data_){data=data_;}
    /// \brief Recorder watching a variable; the stored value is left uninitialised until Save().
    /// \param addr_ address of the variable (not owned; must stay valid)
    recorder(Element * addr_){addr=addr_;}
    ~recorder(){}
    /// \brief Set or change the watched variable.
    /// \param addr_ address of the variable (not owned)
    void    SetAddress(Element * addr_){ addr = addr_; }
    /// \brief The stored value.
    Element GetData(){return data;}
    /// \brief Store the current value of the watched variable unconditionally; does nothing if no address is set.
    void    Save(){ if(!addr){return;} Save(*addr); }
    /// \brief Store the given value unconditionally.
    /// \param data_ value to store
    void    Save(const Element data_){data=data_;}
    /// \brief Store the watched variable if it is larger than the stored value.
    /// \return true if saved
    /// \warning Dereferences the address without checking: crashes if none was set.
    Bool_t  MaxSave();
    /// \brief Same as MaxSave() for a given value, not the watched variable.
    /// \param data_ candidate value
    /// \return true if saved
    Bool_t  MaxSave(const Element data_);
    /// \brief Store the watched variable if it is smaller than the stored value.
    /// \return true if saved
    /// \warning Dereferences the address without checking: crashes if none was set.
    Bool_t  MinSave();
    /// \brief Store \c data_ if it is smaller than the stored value.
    /// \param data_ candidate value
    /// \return true if saved
    Bool_t  MinSave(const Element data_);
    /// \brief Print the watched variable's value and the stored value on \c std::cout (a note is printed if no address is set).
    void    Print();
    protected:
    /// \brief Compare the watched variable with the stored value and save it (shared by MaxSave and MinSave).
    /// \param sgn true for the maximum variant, false for the minimum variant
    /// \return true if saved
    Bool_t   CMPSave(const Bool_t sgn);
    /// \brief Compare \c data_ with the stored value and save it if the test passes.
    /// \param sgn true for the maximum variant, false for the minimum variant
    /// \param data_ candidate value
    /// \return true if saved
    Bool_t   CMPSave(const Bool_t sgn, const Element data_);
    private:
    Element data;               ///< stored value
    Element * addr = nullptr;   ///< watched variable, not owned
};

typedef recorder<Bool_t>   recorderB; ///< recorder of Bool_t
typedef recorder<Int_t>    recorderI; ///< recorder of Int_t
typedef recorder<Long_t>   recorderL; ///< recorder of Long_t
typedef recorder<Float_t>  recorderF; ///< recorder of Float_t
typedef recorder<Double_t> recorderD; ///< recorder of Double_t
#endif