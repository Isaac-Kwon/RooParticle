/// \file runinfo.hpp
/// \brief Run metadata helpers: store how a result file was produced inside the file.

#ifndef H_RUNINFO
#define H_RUNINFO 1

#include "TROOT.h"
#include "TString.h"

/// \brief Run metadata helpers for the test programs.
///
/// The parameters of a run are stored in the current ROOT file (\c gDirectory)
/// as \c TNamed objects (name = key, title = value), so that a result file
/// always tells how it was produced. Read them back with
/// \c file->Get<TNamed>("seed")->GetTitle().
namespace runinfo{
  /// Seed used when no seed is given on the command line (job index 0).
  const UInt_t kDefaultSeed = 65539;

  /// \brief Seed of one job.
  ///
  /// The explicit seed if given (\c seed >= 0), otherwise
  /// \c kDefaultSeed + \c jobIndex, so parallel jobs with different indices
  /// never share a seed and a rerun of the same index is reproducible.
  /// \param jobIndex index of the job (for example the Condor process number)
  /// \param seed explicit seed, or a negative value (default -1) for none; an explicit seed is truncated to 32 bits
  /// \return the seed
  UInt_t SeedFromJob(Int_t jobIndex, Long64_t seed = -1);

  /// \brief Git revision the binary was built from.
  /// \return the short SHA that CMake passes as the macro \c RP_GIT_SHA, or
  /// "unknown" if that macro is not defined (no git, or no \c .git directory)
  TString GitSHA();

  /// \name Write one key/value pair
  /// Writes a \c TNamed(key, value) into the current directory (\c gDirectory)
  /// with the key \c key, overwriting an earlier entry of that name. Does
  /// nothing if there is no current directory. Call before closing the file.
  /// @{

  /// \param key name of the entry
  /// \param value text stored as title
  void Write(const char* key, const char* value);
  /// \param key name of the entry
  /// \param value number, stored as text with 17 significant digits ("%.17g")
  void Write(const char* key, Double_t value);
  /// \param key name of the entry
  /// \param value integer, stored as decimal text
  void Write(const char* key, Long64_t value);
  /// @}

  /// \brief Write the parameters every run shares.
  ///
  /// Writes the keys \c seed, \c dt, \c integrator, \c nevents, \c inspectors
  /// and \c git_sha (the last one is added automatically from GitSHA()) into
  /// the current ROOT directory.
  /// \param seed random seed of the run
  /// \param dt time step in fm/c
  /// \param integrator integrator name (for example event::GetIntegratorName())
  /// \param nevents number of events
  /// \param inspectorCriteria text describing the inspector criteria
  void WriteCommon(UInt_t seed, Double_t dt, const char* integrator,
                   Long64_t nevents, const char* inspectorCriteria);
}

#endif
