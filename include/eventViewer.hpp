/// \file eventViewer.hpp
/// \brief Unfinished, disabled stub for a graphical event viewer.
///
/// \warning Not implemented. Everything in this header is commented out (a
/// GUIparticle with a colour and an eventViewer returning a TCanvas); only the
/// includes remain active, and src/eventViewer.cpp is empty apart from its
/// include. Nothing in the code base uses it.
/// \todo Implement or remove.

#ifndef H_EVENTVIEWER
#define H_EVENTVIEWER 1;

#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "event.hpp"

// class GUIparticle: public particle{
// public:
//   GUIParticle();
//   TColor * GetColor(){return color;}
// protected:
// private:
//   TColor * color = new TColor();
// };

// class eventViewer{
// public:
//   eventViewer(event* event_){tevent = event_;}
//   TCanvas * GetCanvas();
//   event * GetEvent();
// protected:
// private:
//   event * tevent;
// };

#endif
