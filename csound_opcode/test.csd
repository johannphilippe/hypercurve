<CsoundSynthesizer>
<CsOptions>
-o dac
</CsOptions>
<CsInstruments>
sr = 48000
ksmps =  32
nchnls = 2
0dbfs  = 1

gicp = hc_control_point(0.2, 0.6)
gicpp = hc_control_point(0.4, 0.6)

gidiocles = hc_hypercurve(0, 2048, 0, 
                hc_segment(1/2, 1, hc_diocles_curve(1)), 
                hc_segment(1/2, 0, hc_diocles_curve(1)))
gicub = hc_hypercurve(0, 2048, 0, 
                hc_segment(1/2, 1, hc_cubic_curve()), 
                hc_segment(1/2, 0, hc_cubic_curve()))
giquad_bez = hc_hypercurve(0, 2048, 0, 
                hc_segment(1/2, 1, hc_quadratic_bezier_curve( hc_control_point(0.2, 0.6) ) ), 
                hc_segment(1/2, 0, hc_quadratic_bezier_curve( hc_control_point(0.5, 0.2) ) ) )

gicub_bez = hc_hypercurve(0, 2048, 0, 
                hc_segment(1/2, 1, hc_cubic_bezier_curve(hc_control_point(0.2, 0.8), hc_control_point(0.5, 0.3))), 
                hc_segment(1/2, 0, hc_cubic_bezier_curve(hc_control_point(0.8, 0.1), hc_control_point(0.9, 0.7))))

gicm = hc_hypercurve(0, 2048, 0, 
                hc_segment(1/2, 1, hc_catmull_rom_curve(hc_control_point(-1, -5), hc_control_point(2,  4))), 
                hc_segment(1/2, 0, hc_catmull_rom_curve(hc_control_point(-1, -8), hc_control_point(4, 2))))

gipolynomial = hc_hypercurve(0, 2048, 0,
                hc_segment(1, 1, hc_polynomial_curve(5, 3, -1.5, 9)))

gilagrange = hc_hypercurve(0, 2048, 0, 
                hc_segment(1, 1, hc_lagrange_curve(  hc_control_point(0.2, 0.8), hc_control_point(0.4, 0.1), hc_control_point(0.75, 0.5))))


;hc_normalize_y(gipolynomial, 0, 1)
;hc_write_as_png(gipolynomial, "./poly.png", 1, 1)

print(gicub)
print(gidiocles)
giadd = hc_add(gicub, gidiocles)
gisub = hc_sub(gicub_bez, gipolynomial)
gimult = hc_mult(gicub, gidiocles)
gidiv = hc_div(gicm, gipolynomial)

hc_normalize_y(giadd, 0, 1)
hc_normalize_y(gisub, 0, 1)
hc_normalize_y(gimult, 0, 1)
hc_normalize_y(gidiv, 0, 1)


instr 1
        icurve = p4
        kenv = tablei:k( linseg(0, p3, 1), icurve)
        ao = vco2(0.3,  200 + (kenv * 200) ) * kenv
        outs ao,ao
endin

instr 9
        schedule(1, 0, 3, gidiocles)
        schedule(1, 3, 3, gicub)
        schedule(1, 6, 3, giquad_bez)
        schedule(1, 9, 3, gicub_bez)
        schedule(1, 12, 3, gicm)
        schedule(1, 15, 3, gilagrange)

endin

</CsInstruments>
<CsScore>

f 0 z 
i 9 0 0 


</CsScore>
</CsoundSynthesizer>

