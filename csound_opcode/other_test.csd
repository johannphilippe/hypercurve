<CsoundSynthesizer>
<CsOptions>
--opcode-lib=/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug/libcsound_hypercurve.so
-odac
</CsOptions>

<CsInstruments>
sr = 48000
ksmps = 32
nchnls = 2
0dbfs = 1

instr 1
	irnd = jorand()
	print irnd
endin

</CsInstruments>
<CsScore>
i 1 0 5
i 1 + 5
i 1 + 5

</CsScore>
</CsoundSynthesizer>
