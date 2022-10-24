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

segment(frac, y_dest, algo, y_start, rel_x) = environment
{
    abs_diff = abs(y_start - y_dest);
    offset = min(y_start, y_dest);
    scaled = (rel_x) , (1.0 - rel_x) :  select2(y_start > y_dest) : algo : *(abs_diff) : +(offset) ;
    res = scaled;
};

hypercurve(def, y_start_, seg_list, phas) = ba.tabulate(1, proc, def, 0, 1, phas).lin
with {
    seg_count = ba.count(seg_list);

//    fractional_list(indx) = seg_count, gen_frac, indx : rdtable;
//    additional_frac_list(indx) = seg_count, gen_additional, indx : rdtable;

//    gen_frac(size) = ba.take(ba.time, seg_list).frac;
//    gen_additional(size) = sum(n, ba.time, fractional_list); 

    fractional_list = par(n, seg_count, (seg_list : ba.selectn(seg_count,n).frac ) ); 
    real_frac_list = par(n, seg_count,  rint( seg_list : ba.selectn(seg_count, n).frac) );

    additional_frac_list = par(n, seg_count, do_the_sum)
    with {
        do_the_sum(n) = sum(nn, n + 1, fractional_list : ba.selectn(seg_count, nn)); 
    };

    crossed_one(x) = 0, 1 : select2(sum_all > 0) 
    with 
    {
      sum_all = sum(n, seg_count, fc(n));
      fc(n) = (x' < comp) && (x >= comp)
      with { 
        comp = additional_frac_list :  ba.selectn(seg_count, n ); 
      };
    };

        
    proc(x) = which 
    with { 
        cur_seg_idx = ba.counter(crossed_one(x));
        /*
        cur_frac = fractional_list : ba.selectn(seg_count, cur_seg_idx);
        current_seg_factor = 1.0 / cur_frac; 
        to_deduce =  additional_frac_list : ba.selectn(seg_count, cur_seg_idx );
        x_deduced = (x - to_deduce);
        calc_x = x_deduced * current_seg_factor;

        cur_seg  = seg_list : ba.selectn( seg_count,cur_seg_idx + 1);
        y_start = y_start_, old_dest : select2(cur_seg_idx > 0)
        with {
            old_dest = seg_list : ba.selectn(seg_count, cur_seg_idx).y_dest;
        };
        which = y_start, calc_x : cur_seg;
        */
        which = cur_seg_idx;
    };
};


hc = os.phasor(1, 0.5) : hypercurve(4096, 0, (segment(0.5, 1, cubic_curve), segment(0.5, 0, cubic_curve) ));

process = os.osc(150) * 0.1 * hc;
