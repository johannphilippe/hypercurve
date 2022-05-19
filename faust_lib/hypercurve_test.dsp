os = library("oscillators.lib");
hc = library("hypercurve.lib");
si = library("signals.lib");
ve = library("vaeffects.lib");
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
fq_amp = hslider("freq_mult", 0.1, 0, 0.8, 0.001);

res = hc.run(crv, os.phasor(1, 0.3));
process = os.sawtooth(100) : ve.korg35LPF(fq_amp * res + 0.1, 1) : *(amp) : *(res);
