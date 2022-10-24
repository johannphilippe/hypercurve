import("stdfaust.lib");
// MVCE 

diocles_curve(A, x) = process_diocles(x) * COMPENSATION
with 
{
    process_diocles(x) = sqrt( (x*x*x) / (2 * A - x) );
    COMPENSATION = 1.0 / process_diocles(1.0);
};

cubic_curve(x) = process_cubic(x) * COMPENSATION
with 
{
    process_cubic(x) = x*x*x;
    COMPENSATION = 1.0 / process_cubic(1.0);
};

linear(x) = x;

segment(rel_x, y_start, y_dest, algo) =  scaled
with {
    abs_diff = abs(y_start - y_dest);
    offset = min(y_start, y_dest);
    scaled = (rel_x) , (1.0 - rel_x) :  select2(y_start > y_dest) : algo : *(abs_diff) : +(offset) ;
};
//qLog(x) = ba.tabulate(1,log(_),TableSize,MinFreq,MaxFreq,x).cub;
hybrid(def, phas) =  ba.tabulate(1, proc, def, 0.0, 1.0, phas).lin //def, proc(def), int(phas * def) : rdtable
with {
    proc(d) = which
    with {
        x = d ;/// float(def);
        calc_x = (x - to_deduce ) * current_seg_factor;
        to_deduce = 0, 0.5 : select2(x >= 0.5);
        current_seg_factor = 2;
        //current_seg_factor = 1, 2 : select2(x >= 0.5);
        which = dio, cub : select2(x >= 0.5)
        with {
            /*
            dio = segment(calc_x, 0, 1, linear);
            cub = segment(calc_x, 1, 0, linear);
            */
            
            dio = segment(calc_x, 0, 1, diocles_curve(0.51));
            cub = segment(calc_x, 1, 0, cubic_curve);
            /*
            
            dio = 0.51, calc_x : diocles_curve;
            cub = calc_x : cubic_curve;
            */
        };
    };
};

hc = hybrid(16384, os.phasor(1, 1));
bg = hc : vbargraph("hypercurve", 0, 1);
amp = hslider("amp", 1, 0, 1, 0.01) : si.smoo;
process  = os.osc(300 + (500 * hc))  * 0.1 * hc * amp ; //, bg ; 