<CsoundSynthesizer>
<CsOptions>
-odac
</CsOptions>

<CsInstruments>
sr = 48000
ksmps = 32
nchnls = 2
0dbfs = 1

instr 1
	puts("diocles", 1)
	idio = diocles_curve(1)
	puts("seg", 1)
	iseg1 = segment(1, 1, idio)
	puts("crv", 1)
	icrv = hypercurve(2048, 0, iseg1)
	kcrv = run_hypercurve(icrv, phasor:k(0.5))
	ao = vco2(0.3, 400) * kcrv
	outs ao,ao
endin

</CsInstruments>
<CsScore>
i 1 0 5
i 1 + 5
i 1 + 5

</CsScore>
</CsoundSynthesizer>
