function g = design_chirplet(f,sp,signal)


clear g
g.center_frequency = f; %phase_f
g.fractional_bandwidth = 0.25;
g.chirp_rate = 0;
g = make_chirplet('chirplet_structure', g, 'signal_parameters', sp);