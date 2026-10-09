//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/ChannelBlend.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Registry/Registry.h"

namespace genie {
namespace utils {
namespace channelblend {

namespace {

  const char * kKeyEmin = "ChannelBlend-Emin";
  const char * kKeyEmax = "ChannelBlend-Emax";

  // The window is read from the tune's global parameter list the first time
  // that it is needed. The tune must already have been built at that point
  // (it has been, by the time that a driver sums cross sections).
  struct Window {
    bool   active = false;
    double emin   = 0.;
    double emax   = 0.;
  };

  const Window & GetWindow(void)
  {
    static bool   loaded = false;
    static Window win;
    if (loaded) return win;

    const Registry * gc = AlgConfigPool::Instance()->GlobalParameterList();
    if (gc && gc->Exists(kKeyEmin) && gc->Exists(kKeyEmax)) {
      win.emin = gc->GetDouble(kKeyEmin);
      win.emax = gc->GetDouble(kKeyEmax);
      win.active = (win.emax > win.emin);
      if (!win.active) {
        LOG("ChannelBlend", pERROR)
          << "Ignoring the MARLEY blend window: need " << kKeyEmax << " > "
          << kKeyEmin << ", got [" << win.emin << ", " << win.emax << "] GeV";
      } else {
        LOG("ChannelBlend", pNOTICE)
          << "MARLEY blend is active: MARLEY below " << win.emin
          << " GeV, GENIE above " << win.emax << " GeV, linear in between";
      }
    }
    loaded = true;
    return win;
  }

} // anonymous namespace

//____________________________________________________________________________
bool   IsActive(void) { return GetWindow().active; }
double Emin    (void) { return GetWindow().emin;   }
double Emax    (void) { return GetWindow().emax;   }
//____________________________________________________________________________
double MarleyProbability(double E)
{
  const Window & w = GetWindow();
  if (!w.active) return 0.;
  if (E <= w.emin) return 1.;
  if (E >= w.emax) return 0.;
  return (w.emax - E) / (w.emax - w.emin);
}
//____________________________________________________________________________
double Weight(const Interaction & in, const InteractionList & ilst, double E)
{
  if (!IsActive()) return 1.;

  const double pM = MarleyProbability(E);
  const ProcessInfo & pi = in.ProcInfo();

  if (pi.IsMarley()) return pM;

  // Scattering off atomic electrons (nu-e elastic, IMD, ...) and the Glashow
  // resonance are not part of the MARLEY-CC / MARLEY-NC channels (MARLEY
  // simulates reactions on the nucleus), so they are never suppressed.
  if (pi.IsElectronScattering() || pi.IsGlashowResonance()) return 1.;

  // A native GENIE channel is only suppressed if MARLEY covers the same type
  // of interaction (CC, NC, ...) for this initial state.
  const InteractionType_t itype = pi.InteractionTypeId();
  for (InteractionList::const_iterator it = ilst.begin(); it != ilst.end(); ++it) {
    const ProcessInfo & other = (*it)->ProcInfo();
    if (other.IsMarley() && other.InteractionTypeId() == itype) return 1. - pM;
  }
  return 1.;
}
//____________________________________________________________________________

} // channelblend namespace
} // utils namespace
} // genie namespace
