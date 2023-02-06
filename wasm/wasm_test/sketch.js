/*jshint esversion: 11 */

const CSOUND_SOURCE = "https://unpkg.com/@csound/browser@6.18.6/dist/csound.js";
const HYPERCURVE_SOURCE = "hypercurve.wasm"


const start = async () => {  
  const { Csound } = await import(CSOUND_SOURCE);
  const csound = await Csound({ 
    useWorker: false, 
    useSAB: false,
    withPlugins: [HYPERCURVE_SOURCE]
  });
  await csound.setOption("-odac");
  await csound.compileCsdText(`
</CsOptions>
    -odac
<CsInstruments>

sr = 44100
ksmps = 32
nchnls = 2
0dbfs  = 1

gicrv = hc_hypercurve(0, 2048, 0, hc_segment(0.5, 1, hc_cubic_curve()), hc_segment(0.5, 0, hc_cubic_curve()))

instr 1
  ;kcrv = linseg(0,p3/2, 1, p3/2, 0)
  kcrv = tablei:k(linseg:k(0,p3,1), gicrv, 1)
  ao = oscili(0.3, 300) * kcrv
  outs ao,ao
endin

</CsInstruments>
<CsScore>
i1 0 10
e
</CsScore>
</CsoundSynthesizer>
              `);
  await csound.start();
  
}

const triggerEvent = "ontouchstart" in document.documentElement ? 
      "touchend" : "click";

document.querySelector("#play_test").addEventListener(triggerEvent, start);

