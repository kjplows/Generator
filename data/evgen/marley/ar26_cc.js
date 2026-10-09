// MARLEY job file for the MARLEY-CC channel of the GENIE AR26 tune
// (AR26 = AR25 + MARLEY). Used through genie::MarleyInterface/AR26-CC.
//
// Only 40Ar nu_e CC reactions are listed, so that this MARLEY generator can only
// produce CC events. The total cross section that GENIE sees for the MARLEY-CC
// channel is the total of the reactions below.
//
// - Continuum transitions: HF-CRPA
// - Discrete transitions:  Bhattacharya et al. (2009)
// This is the recommended MARLEY v2 combination.
//
// The neutrino source below is a placeholder. GENIE supplies the probe
// energy and direction for every event. (MARLEY only requires it so that at
// least one reaction can be matched to a source neutrino at initialization.)
{
  target: {
    nuclides: [ 1000180400 ],
    atom_fractions: [ 1.0 ],
  },

  reactions: [
    "ve40ArCC_HF-CRPA.react",
    "ve40ArCC_Bhattacharya2009-Discrete.react",
  ],

  source: {
    type: "dar",
    neutrino: "ve",
  },

  direction: { x: 0.0, y: 0.0, z: 1.0 },

  // MARLEY simulates the de-excitation of the residual nucleus (HF decays and
  // tabulated gamma cascades) as part of create_event()
  do_deexcitations: true,
}
