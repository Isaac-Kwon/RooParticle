/// \file path.hpp
/// \brief Recorded trajectory of one particle (\ref path) and distance-of-closest-approach helpers.

#ifndef H_PATH
#define H_PATH 1;

#include "TROOT.h"
#include "RPConfig.hpp"
#include "TVector3.h"
#include "iostream"

#include "list"

using namespace std;

/// \brief Time-ordered list of (position, velocity, time) samples of one particle.
///
/// A \ref particle owns one path when its recording is switched on; every
/// recorded step appends one record with the particle's position (fm),
/// velocity (c) and its own clock \c etime (fm/c). The three internal vectors
/// always have the same length, record \c i is the \c i-th stored sample (not
/// the \c i-th simulation step: see particle::SetRecordPeriod).
///
/// \warning With a record period n > 1 only every n-th step is stored, so
/// GetLast() / GetLastX() is NOT the final state of the particle (use the
/// particle itself for that). The DCA search also sees only the stored
/// samples, so a coarse record period makes the discrete minimum coarse
/// (the interpolation of GetDCA() only partly compensates).
///
/// \note An invincible (fixed) particle does not advance its own clock, so the
/// times in its path stay constant; such a path does not share a time grid
/// with the path of a moving particle (see the path-vs-path GetDCA()). For a
/// fixed centre use the point overload instead.
///
/// \note The getters never throw. For an out-of-range index or an empty path
/// they print a message on \c std::cerr and return a null vector / 0.
class path{
public:
  /// \brief Create an empty path.
  path();
  /// \brief Copy constructor: copies all records (deep copy of the three vectors).
  path(const path& path_);
  /// \brief Destructor: clears the records.
  ~path();
  /// \brief Append one record at the end.
  /// \param x_ position in fm
  /// \param v_ velocity in c
  /// \param time_ time in fm/c (the particle's own clock)
  /// \note The definition in path.cpp declares a default time of -1, but the
  /// header has none, so callers must always pass the time.
  void AddRecord(TVector3 x_, TVector3 v_, Double_t time_);

  /// \brief Number of stored records (size of the position vector).
  Int_t    GetRecordSize(){return pos.size();}

  /// \brief Position of record \c i (fm).
  /// \param i record index, 0 .. size-1
  /// \return the position, or a null vector (and a message on \c std::cerr)
  /// if the path is empty or \c i is out of range (a negative \c i also counts)
  TVector3 GetRecordX(Int_t i);
  /// \brief Velocity of record \c i (c); same failure behaviour as GetRecordX().
  TVector3 GetRecordV(Int_t i);
  /// \brief Time of record \c i (fm/c).
  /// \return the time, or 0 (and a message on \c std::cerr) if the path is
  /// empty or \c i is out of range. Note that 0 is also a valid time.
  Double_t GetRecordT(Int_t i);
  /// \brief Position, velocity and time of record \c i in one call.
  /// \param i record index
  /// \param[out] x_ position (fm)
  /// \param[out] v_ velocity (c)
  /// \param[out] t_ time (fm/c)
  /// Each output falls back individually as in GetRecordX(), GetRecordV(), GetRecordT().
  void     GetRecord(Int_t i, TVector3 & x_, TVector3 & v_, Double_t & t_);

  /// \brief Position of the first record (fm); null vector plus message if the path is empty.
  TVector3 GetInitX();
  /// \brief Velocity of the first record (c); null vector plus message if the path is empty.
  TVector3 GetInitV();
  /// \brief Time of the first record (fm/c); 0 plus message if the path is empty.
  Double_t GetInitT();
  /// \brief First record as position (fm), velocity (c) and time (fm/c).
  /// \param[out] x_ position
  /// \param[out] v_ velocity
  /// \param[out] t_ time
  void     GetInit(TVector3 & x_, TVector3 & v_, Double_t & t_);

  /// \brief Position of the last stored record (fm); null vector plus message if empty.
  /// \warning Not the final state of the particle when the record period is > 1.
  TVector3 GetLastX();
  /// \brief Velocity of the last stored record (c); null vector plus message if empty.
  TVector3 GetLastV();
  /// \brief Time of the last stored record (fm/c); 0 plus message if empty.
  Double_t GetLastT();
  /// \brief Last stored record as position (fm), velocity (c) and time (fm/c).
  /// \param[out] x_ position
  /// \param[out] v_ velocity
  /// \param[out] t_ time
  void     GetLast(TVector3 & x_, TVector3 & v_, Double_t & t_);

  /// \brief Distance of closest approach of this path to a fixed point.
  ///
  /// The distance squared \f$d_i^2 = |\vec x_i - \vec p|^2\f$ is evaluated at
  /// every stored record and the smallest one (first one on a tie) is taken.
  /// With \c interpolate the result is refined: a parabola through the
  /// distance squared of the three records around the discrete minimum
  /// (indices k-1, k, k+1; the time grid need not be uniform) gives the
  /// sub-step minimum and the time at which it occurs. The refined value is
  /// \f$\sqrt{d^2_{\min}}\f$ with \f$d^2_{\min}\f$ the parabola's vertex value.
  /// The discrete values are kept instead if the minimum is the first or last
  /// record, if the times around it are not strictly increasing, if the
  /// parabola does not open upwards, or if its vertex lies outside
  /// [t(k-1), t(k+1)]. Without \c interpolate the discrete minimum is returned.
  ///
  /// \param pos_ the point (fm); only the path moves, the point is fixed
  /// \param[out] DCA distance of closest approach (fm), or -1 on failure
  /// \param[out] outtime time of the closest approach (fm/c, the path's own
  ///        clock), or -1 on failure
  /// \param verbose if true, print the index, DCA and time on \c std::cerr
  /// \param interpolate if true (default), refine by the parabola
  ///
  /// \note Failure: if the path is empty (or its position and time vectors
  /// differ in length), both outputs are set to the sentinel -1 and a message
  /// goes to \c std::cerr. The sentinel -1 is distinguishable from a result
  /// because a real distance is never negative.
  /// \note The path only knows its stored samples. With a record period > 1
  /// the true minimum may be missed between samples.
  void     GetDCA(const TVector3 pos_, Double_t & DCA, Double_t & outtime, Bool_t verbose=false, Bool_t interpolate=true);
  /// \brief Distance of closest approach between this path and another path.
  ///
  /// The two paths are compared record by record, so they must lie on the
  /// same time grid: the same number of records and equal times in every
  /// record (exact floating-point equality). The relative distance squared
  /// \f$|\vec x^{(a)}_i - \vec x^{(b)}_i|^2\f$ is minimised exactly as in
  /// GetDCA(const TVector3, Double_t&, Double_t&, Bool_t, Bool_t) and refined
  /// by the same parabola when \c interpolate is true. The time is taken from
  /// this path.
  ///
  /// \param timepath the other path, not owned; must not be null
  /// \param[out] DCA distance of closest approach (fm), or -1 on failure
  /// \param[out] outtime time of the closest approach (fm/c), or -1 on failure
  /// \param verbose if true, print the index, DCA and time on \c std::cerr
  /// \param interpolate if true (default), refine by the parabola
  ///
  /// \note Failure (both outputs set to the sentinel -1, message on
  /// \c std::cerr): \c timepath is null; the paths differ in length or in any
  /// time; this path is empty.
  /// \warning Two paths from particles with different clocks (e.g. an
  /// invincible particle, whose clock does not advance) or with different
  /// record periods do not share a time grid and give -1.
  void     GetDCA(path* timepath, Double_t & DCA, Double_t & outtime, Bool_t verbose=false, Bool_t interpolate=true);
  /// \brief Number of stored time stamps (as a 64-bit integer); equals GetRecordSize().
  Long64_t GetMaxNumber(){return t.size();}
  // Double_t GetDCA(const TVector3 pos);
  // Double_t GetDCA(path* timepath);
protected:
  /// \brief True if \c target has the same number of records; otherwise prints a message on \c std::cerr and returns false.
  Bool_t lengthChk(path* target);
  /// \brief True if \c target has the same length and exactly equal times in every record; prints a message on \c std::cerr and returns false otherwise.
  Bool_t timeChk(path* target);

private:
  // Int_t writingCursor = 0;
  std::vector<TVector3>  pos; ///< positions (fm)
  std::vector<TVector3>  vel; ///< velocities (c)
  std::vector<Double_t>  t;   ///< times (fm/c)
};

#endif
