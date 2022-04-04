<CsoundSynthesizer>
<CsOptions>
-o dac
--opcode-lib=./csound_opcode/libcsound_hypercurve.so
</CsOptions>
<CsInstruments>
sr = 48000
ksmps =  32
nchnls = 2
0dbfs  = 1

gicp = control_point(0.2, 0.6)
gicpp = control_point(0.4, 0.6)

gidiocles = hypercurve(2048, 0, 
		segment(1/2, 1, diocles_curve(1)), 
		segment(1/2, 0, diocles_curve(1)))
gicub = hypercurve(2048, 0, 
		segment(1/2, 1, cubic_curve()), 
		segment(1/2, 0, cubic_curve()))
giquad_bez = hypercurve(2048, 0, 
		segment(1/2, 1, quadratic_bezier_curve( control_point(0.2, 0.6) ) ), 
		segment(1/2, 0, quadratic_bezier_curve( control_point(0.5, 0.2) ) ) )

gicub_bez = hypercurve(2048, 0, 
		segment(1/2, 1, cubic_bezier_curve(control_point(0.2, 0.8), control_point(0.5, 0.3))), 
		segment(1/2, 0, cubic_bezier_curve(control_point(0.8, 0.1), control_point(0.9, 0.7))))

gicm = hypercurve(2048, 0, 
		segment(1/2, 1, catmull_rom_curve(control_point(-1, -5), control_point(2,  4))), 
		segment(1/2, 0, catmull_rom_curve(control_point(-1, -8), control_point(4, 2))))
instr 1
	icurve = p4
	kenv = run_hypercurve(icurve, linseg(0, p3, 1))

	ao = vco2(0.3,  200 + (kenv * 200) ) * kenv
	outs ao,ao
endin

instr 9
	schedule(1, 0, 3, gidiocles)
	schedule(1, 3, 3, gicub)
	schedule(1, 6, 3, giquad_bez)
	schedule(1, 9, 3, gicub_bez)
	schedule(1, 12, 3, gicm)

endin

</CsInstruments>
<CsScore>

f 0 z 
i 9 0 0 


</CsScore>
</CsoundSynthesizer>

