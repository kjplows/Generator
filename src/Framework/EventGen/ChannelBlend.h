//____________________________________________________________________________
/*!

\namespace  genie::utils::channelblend

\brief      Energy-dependent weights used to blend a "native GENIE" description
            of a given interaction type with a MARLEY-based one in the region
            where both are valid.

\details    The blend is defined at the level of the cross sections of the
            channels that the event generation driver sums over (see
            GEVGDriver::XSecSum() and PhysInteractionSelector).

            Let p_M(E) be the MARLEY probability

                p_M(E) = 1                          E <= Emin
                p_M(E) = (Emax - E)/(Emax - Emin)   Emin < E < Emax
                p_M(E) = 0                          E >= Emax

            Channels with ProcessInfo::IsMarley() are weighted by p_M(E).
            Every other channel is weighted by 1 - p_M(E), *provided that* the
            interaction list also contains a MARLEY channel with the same
            interaction type (CC, NC, ...) for the same initial state.
            Otherwise the weight is 1 and nothing changes (so, for example, a
            MARLEY-CC channel only modifies the CC cross sections).
            Scattering off atomic electrons and the Glashow resonance are
            never suppressed, since MARLEY-CC/NC describe reactions on the
            nucleus only.

            Note that this is equivalent to deciding which model describes
            the neutrino *before* the interaction is accepted (MARLEY with
            probability p_M), and then letting the chosen model's own cross
            section accept or reject the neutrino.

            Because the weights multiply the cross sections themselves, the
            total cross section seen by the flux-driven neutrino acceptance
            (GMCJDriver) is the blended one,

                sigma_blend = p_M sigma_MARLEY + (1 - p_M) sum(sigma_GENIE),

            and the MARLEY channel is subsequently selected with probability
            p_M sigma_MARLEY / sigma_blend by the interaction selector.

            The blend is only active if the global parameter list of the tune
            defines "ChannelBlend-Emin" and "ChannelBlend-Emax" (GeV).
            Otherwise all weights are 1, so tunes without a MARLEY channel
            are unaffected.

\author     MARLEY+GENIE hybrid (AR26)

\created    October 2026
*/
//____________________________________________________________________________

#ifndef _CHANNEL_BLEND_H_
#define _CHANNEL_BLEND_H_

namespace genie {

class Interaction;
class InteractionList;

namespace utils {
namespace channelblend {

  /// True if the tune defines the blend window
  bool IsActive(void);

  /// Lower edge of the blend window (GeV), pure MARLEY below it
  double Emin(void);

  /// Upper edge of the blend window (GeV), pure GENIE above it
  double Emax(void);

  /// Probability p_M(E) that a neutrino of energy E (GeV) is handled by MARLEY
  double MarleyProbability(double E);

  /// Weight to apply to the cross section of the interaction `in` at probe
  /// energy E (GeV). `ilst` is the full interaction list of the driver.
  double Weight(const Interaction & in, const InteractionList & ilst, double E);

} // channelblend namespace
} // utils namespace
} // genie namespace

#endif // _CHANNEL_BLEND_H_
