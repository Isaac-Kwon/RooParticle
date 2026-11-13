/// \file eventVolume.hpp
/// \brief Spatial region attached to an event (\ref eventVolume, \ref eventVGeneral).

#ifndef H_EVENTVOLUME
#define H_EVENTVOLUME 1;

#include "TROOT.h"
#include "TVector3.h"

/// \brief Base class of event volumes: the region in which an event is "alive".
///
/// The base class is the whole space (Inside() is always true). The volume
/// never affects the dynamics. It matters only through the inspectors
/// \c "VOUT" and \c "VOUTA" (event::IsInsideVolume), which ask whether the
/// non-invincible particles have left it, and an event without a volume
/// regards everything as inside. The event does not own the volume
/// (event::SetVolume): create it with \c new or on the stack and keep it alive
/// as long as the event and its copies use it.
class eventVolume{
public:
  /// \brief Create the volume (nothing to initialise).
  eventVolume();
  /// \brief Is the point inside the volume?
  /// \param position point in fm
  /// \return true (the base class contains every point); derived classes override
  virtual Bool_t Inside(TVector3 position){return true;}
protected:
private:
};

/// \brief Axis-aligned box.
///
/// The box is the closed set \f$x_{min}\le x\le x_{max}\f$ and likewise in y
/// and z (fm): a point exactly on a face counts as inside. The constructors
/// do not check that min is smaller than max; a reversed or empty range
/// makes Inside() false for every point.
class eventVGeneral: public eventVolume{
public:
  /// \brief Box from its two corners.
  /// \param xmin lower x bound (fm)
  /// \param xmax upper x bound (fm)
  /// \param ymin lower y bound (fm)
  /// \param ymax upper y bound (fm)
  /// \param zmin lower z bound (fm)
  /// \param zmax upper z bound (fm)
  eventVGeneral(Double_t xmin, Double_t xmax,
                Double_t ymin, Double_t ymax,
                Double_t zmin, Double_t zmax);
  /// \brief Box centred at the origin with the given edge lengths.
  ///
  /// Equivalent to the six-bound constructor with bounds \f$\pm\mathrm{len}/2\f$
  /// in each axis, so eventVGeneral(40.,40.,40.) is the cube
  /// \f$|x|,|y|,|z|\le 20\f$ fm.
  /// \param xlen full edge length along x (fm)
  /// \param ylen full edge length along y (fm)
  /// \param zlen full edge length along z (fm)
  eventVGeneral(Double_t xlen, Double_t ylen, Double_t zlen);
  /// \brief Is the point inside the box (boundary included)?
  /// \param position point in fm
  /// \return true if all three coordinates lie within the bounds, boundary
  /// included; false otherwise
  Bool_t Inside(TVector3 position);
protected:
private:
  TVector3 bottom; ///< lower corner (xmin, ymin, zmin) in fm
  TVector3 top;    ///< upper corner (xmax, ymax, zmax) in fm
};

#endif
