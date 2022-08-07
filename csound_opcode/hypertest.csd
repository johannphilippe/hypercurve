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
	icrv = hc_hypercurve(0, 16384, 0, 
			hc_segment(1/8, 1, hc_gauss_curve(10, 0.5)), 
			hc_segment(1/8, 0.5, hc_blackman_curve()),
			hc_segment(4/8, 0.2, hc_catmull_rom_curve( hc_control_point(-2, -1), hc_control_point(2, 2))),
			hc_segment(2/8, 0, hc_cissoid_curve(1)))


	kenv = tablei:k(linseg:k(0, p3, 1), icrv, 1)
	ao = lowresx(vco2(0.3, p4), p4 + (kenv * p4), (kenv * 0.3) + 0.1) * kenv

	outs ao, ao
endin

gicrv = hc_hypercurve(0, 16384, 0, 
		hc_segment(1/8, 1, hc_cissoid_curve(2)), 
		hc_segment(1/8, 0.5, hc_power_curve(5)),
		hc_segment(4/8, 0.2, hc_cubic_bezier_curve( hc_control_point(0.3, 0.8), hc_control_point(0.9, 0.1))),
		hc_segment(2/8, 0, hc_hanning_curve()))

giwav = hc_hypercurve(0, 16384, 0, 
		hc_segment(1/3, 1, hc_cissoid_curve(1)), 
		hc_segment(1/3, -1, hc_hamming_curve()), 
		hc_segment(1/3, 0, hc_quadratic_bezier_curve(hc_control_point(0.3, 0.8))))

instr 2
	icrv = gicrv
	iwav = giwav
	kenv = tablei:k(linseg:k(0, p3, 1), icrv, 1)
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

