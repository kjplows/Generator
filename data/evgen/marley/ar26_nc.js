// MARLEY job file for the MARLEY-NC channel of the GENIE AR26 tune
// (AR26 = AR25 + MARLEY). Used through genie::MarleyInterface/AR26-NC.
//
// Only coherent elastic neutrino-nucleus scattering (CEvNS) on 40Ar is listed,
// so that this MARLEY generator can only produce NC events of that type. The
// reaction is registered by MARLEY for all four (anti)neutrino flavours.
// (Neutrino-electron elastic scattering is deliberately NOT included here: GENIE
// keeps simulating it natively with its NUE-EL channel.)
//
// The neutrino source below is a placeholder. GENIE supplies the probe
// energy and direction for every event.
{
  target: {
    nuclides: [ 1000180400 ],
    atom_fractions: [ 1.0 ],
  },

  reactions: [
    "CEvNS40Ar.react",
  ],

  source: {
    type: "dar",
    neutrino: "ve",
  },

  direction: { x: 0.0, y: 0.0, z: 1.0 },

  do_deexcitations: true,
}
