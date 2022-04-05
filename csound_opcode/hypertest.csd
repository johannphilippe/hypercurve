<CsoundSynthesizer>
<CsOptions>
-o dac
--opcode-lib=./libcsound_hypercurve.so
</CsOptions>
<CsInstruments>
sr = 48000
ksmps =  32
nchnls = 2
0dbfs  = 1



instr 1
	icrv = hypercurve(16384, 0, 
			segment(1/8, 1, gauss_curve(10, 0.5)), 
			segment(1/8, 0.5, blackman_curve()),
			segment(4/8, 0.2, catmull_rom_curve( control_point(-2, -1), control_point(2, 2))),
			segment(2/8, 0, cissoid_curve(1)))

	kenv = run_hypercurve(icrv, linseg:k(0, p3, 1))
	ao = lowresx(vco2(0.3, p4), p4 + (kenv * p4), (kenv * 0.7) + 0.1) * kenv

	outs ao, ao
endin

gicrv = hypercurve(16384, 0, 
		segment(1/8, 1, cissoid_curve(2)), 
		segment(1/8, 0.5, power_curve(5)),
		segment(4/8, 0.2, cubic_bezier_curve( control_point(0.3, 0.8), control_point(0.9, 0.1))),
		segment(2/8, 0, hanning_curve()))

giwav = hypercurve(16384, 0, 
		segment(1/3, 1, cissoid_curve(1)), 
		segment(1/3, -1, hamming_curve()), 
		segment(1/3, 0, quadratic_bezier_curve(control_point(0.3, 0.8))))

instr 2
	icrv = gicrv
	iwav = giwav
	kenv = run_hypercurve(icrv, linseg:k(0, p3, 1))
	//	ares = run_hypercurve(iwav, phasor:a(p4))
	ares = vco2(1, p4, 0) 
	ares *= kenv * 0.3
	outs ares, ares
endin

</CsInstruments>
<CsScore>

i 1 0 2 100
i 1 0 2 200

i 1 1.5 3 120
i 1 1.5 3 160

i 1 3 8 150
i 1 3 8 220
i 1 6 10 180
i 1 6 10 200
i 1 6 10 270

i 1 14 2 300
i 1 14 2 250
i 1 14.5 2 500
i 1 14.5 2 100
i 1 16 2 260

i 2 0 10 80

</CsScore>
</CsoundSynthesizer>

