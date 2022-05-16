os = library("oscillators.lib");
hc = library("hypercurve.lib");
crv = hc.hypercurve(2048, 0, (hc.segment(0.3, 1, hc.diocles_curve(1)), (hc.segment(1/3, 1, hc.linear_curve), hc.segment(0.3, 0, hc.diocles_curve(0.6)))));
res = hc.run(crv, os.phasor(1, 1));
process = os.sawtooth(440) * res;
