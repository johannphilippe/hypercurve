os = library("oscillators.lib");
hc = library("hypercurve.lib");
si = library("signals.lib");
ve = library("vaeffects.lib");
sp = library("spats.lib");
no = library("noises.lib");
ba = library("basics.lib");
// Simple example
// crv = hc.hypercurve(2048, 0, (hc.segment(0.3, 1, hc.diocles_curve(1)), (hc.segment(1/3, 1, hc.linear_curve), hc.segment(0.3, 0, hc.diocles_curve(0.6)))));

// Spline example
crv = hc.hypercurve(2048, 0,
     (hc.segment(1/2, 1, hc.cubic_bezier_curve(hc.control_point(0.45, 0.02), hc.control_point(0.55, 0.98))),
      hc.segment(1/2, 0, hc.cubic_curve)
      ));

//crv = hc.hypercurve(2048, 0,
//     (hc.segment(1/2, 1, hc.quadratic_bezier_curve(hc.control_point(0.45, 0.02)))));

//crv = hc.hypercurve(2048, 0,
//     (hc.segment(7/8, 1, hc.catmull_rom_spline_curve(hc.control_point(-1, -4), hc.control_point(2, 8))),
//      hc.segment(1/8, 0, hc.gaussian_curve(1.2, 0.5)))
//     );


//crv = hc.hypercurve(2048, 0,
//     (hc.segment(3/4, 1, hc.cubic_spline_curve( (
//                hc.control_point(0.25, 0.1),
//                hc.control_point(0.45, 0.15),
//                hc.control_point(0.55, 0.8),
//                hc.control_point(0.90, 0.1))
//                )),
//      hc.segment(1/4, 0, hc.catenary_curve(0.1))
//              ));
// Sometimes, cubic spline needs to be rescaled
// hc.normalize_y(crv, 0, 1)

//crv = hc.hypercurve(4096, 0, (hc.segment(1, 1, hc.linear_curve)));

//crv = hc.hypercurve(4096, 0, ( hc.segment(1, 1, hc.lagrange_polynomial_curve( ( hc.control_point(0.2, 0.8), hc.control_point(0.4, 0.1))))));

amp = hslider("amp", 0.3, 0, 1, 0.01) : si.smoo;
fq_amp = hslider("freq_mult", 0.1, 0, 0.7, 0.001);
phas_mult = hslider("speed", 0.1, 0.01, 1, 0.01);

snt(freq, ph_freq, PHASE, amp_mult, index , sigL, sigR) = os.sawtooth(freq) : ve.korg35LPF(fq_amp * env + 0.2, 2) : *(amp) : *(env) : *(amp_mult) : sp.panner(pan) : +(sigL), +(sigR)
with {
      env = hc.run(crv, os.hsp_phasor(1,ph_freq * phas_mult, 0, PHASE));
      pan = 1-env, 0+env : select2(index % 2 == 0);
};

N_OSC = 16;

rnd(n) = no.noise + 0.01 : ba.sAndH(os.impulse) : abs;
process = 0, 0 : seq(n, N_OSC, snt( 100 + (500 * rnd(n) * n), (n+1) / N_OSC, (n+1) / N_OSC, 1 / N_OSC, n ));
