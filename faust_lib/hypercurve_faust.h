/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef HYPERCURVE_FAUST_H
#define HYPERCURVE_FAUST_H

#include"../src/hypercurve.h"
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
// Glue code for Faust hypercuve.
////////////////////////////////////////////////////////////////////////////////


using namespace hypercurve;
extern "C"
{
static increment_map< std::shared_ptr<curve_base> > faust_curve_base_map;
static increment_map< std::shared_ptr<control_point> > faust_control_point_map;
static increment_map< std::shared_ptr<segment> > faust_segment_map;
static increment_map< std::shared_ptr<curve> > faust_curve_map;

// Operations on curves
static int hc_scale(int crv, double min, double max)
{
    faust_curve_map[crv]->scale(min, max);
    return crv;
}

// Just inverts a curve base on its linear axis (symetry on linear axis)
static int hc_invert_curve_base(int cb)
{
    vinvert(faust_curve_base_map[cb]);
    return cb;
}

static int hc_mirror_curve_base(int cb)
{
    mirror(faust_curve_base_map[cb]);
    return cb;
}

// Operators on hypercurves
static int hc_add(int h1, int h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] + *faust_curve_map[h2]));
}
static int hc_sub(int h1, int h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] - *faust_curve_map[h2]));
}
static int hc_mult(int h1, int h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] * *faust_curve_map[h2]));
}
static int hc_div(int h1, int h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] / *faust_curve_map[h2]));
}
static int hc_addn(int h1, double h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] + h2));
}
static int hc_subn(int h1, double h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] - h2));
}
static int hc_multn(int h1, double h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] * h2));
}
static int hc_divn(int h1, double h2)
{
    return faust_curve_map.map(share( *faust_curve_map[h1] / h2));
}
// Curve bases functions

static const int max_segments = 64;

static int hc_diocles_curve(double a) {return faust_curve_base_map.map(share(diocles_curve(a)));}
static int hc_linear_curve(int) {return faust_curve_base_map.map(share(linear_curve()));}
static int hc_cubic_curve(int) {return faust_curve_base_map.map(share(cubic_curve()));}
static int hc_logarithmic_curve(int) {return faust_curve_base_map.map(share(logarithmic_curve()));}
static int hc_exponential_curve(int) {return faust_curve_base_map.map(share(exponential_curve()));}
static int hc_power_curve(double exponent) {return faust_curve_base_map.map(share(power_curve(exponent)));}
static int hc_hanning_curve(int) {return faust_curve_base_map.map(share(hanning_curve()));}
static int hc_hamming_curve(int) {return faust_curve_base_map.map(share(hamming_curve()));}
static int hc_blackman_curve(int) {return faust_curve_base_map.map(share(blackman_curve()));}
static int hc_gaussian_curve(double A, double c) {return faust_curve_base_map.map(share(gaussian_curve(A,c)));}
static int hc_toxoid_curve(double a) {return faust_curve_base_map.map(share(toxoid_curve(a)));}
static int hc_catenary_curve(double a) {return faust_curve_base_map.map(share(catenary_curve(a)));}
static int hc_tightrope_walker_curve(double a, double b) {return faust_curve_base_map.map(share(tightrope_walker_curve(a, b)));}

static int hc_typed_curve(double t) {return faust_curve_base_map.map(share(typed_curve(t)));}
static int hc_mouth_curve(int) {return faust_curve_base_map.map(share(mouth_curve()));}
static int hc_bicorn_curve(int boolean) {return faust_curve_base_map.map(share(bicorn_curve(boolean != 0)));}

static int hc_quadratic_bezier_curve(int cp) {
    return faust_curve_base_map.map(share(quadratic_bezier_curve(*faust_control_point_map[cp])));
}
static int hc_cubic_bezier_curve(int cp1, int cp2) {
    return faust_curve_base_map.map(share(cubic_bezier_curve(*faust_control_point_map[cp1], *faust_control_point_map[cp2])));
}
static int hc_catmull_rom_spline_curve(int cp1, int cp2) {
    return faust_curve_base_map.map(share(catmull_rom_spline_curve(*faust_control_point_map[cp1], *faust_control_point_map[cp2])));
}

static int hc_cubic_spline_curve_varg(int n_args, ...)
{
    int i, arg;
    va_list ap;

    std::vector< control_point > cps;
    va_start(ap, n_args);
    for(i = 0; i < n_args; i++) {
        arg = va_arg(ap, int);
        if(arg <= 0) continue;
        cps.push_back(*faust_control_point_map[arg]) ;
    }
    va_end(ap);
    return faust_curve_base_map.map(share(cubic_spline_curve(cps)));
}

static int hc_lagrange_curve_varg(int n_args, ...)
{
    int i, arg;
    va_list ap;

    std::vector< control_point > cps;
    va_start(ap, n_args);
    for(i = 0; i < n_args; i++) {
        arg = va_arg(ap, int);
        if(arg <= 0) continue;
        cps.push_back(*faust_control_point_map[arg]) ;
    }
    va_end(ap);
    return faust_curve_base_map.map(share(lagrange_polynomial_curve(cps)));
}


static int hc_polynomial_curve_varg(int n_args, ...)
{
    int i, arg;
    va_list ap;

    std::vector< double > cps;
    va_start(ap, n_args);
    for(i = 0; i < n_args; i++) {
        arg = va_arg(ap, double);
        if(arg > 999) continue;
        cps.push_back(arg) ;
    }
    va_end(ap);
    return faust_curve_base_map.map(share(polynomial_curve(cps)));
}

static int hc_cubic_spline_curve(
                    int s1 = -1, int s2 = -1, int s3 = -1, int s4 = -1,
                    int s5 = -1, int s6 = -1, int s7 = -1, int s8 = -1,
                    int s9 = -1, int s10 = -1, int s11 = -1, int s12 = -1,
                    int s13 = -1, int s14 = -1, int s15 = -1, int s16 = -1,
                    int s17 = -1, int s18 = -1, int s19 = -1, int s20 = -1,
                    int s21 = -1, int s22 = -1, int s23 = -1, int s24 = -1,
                    int s25 = -1, int s26 = -1, int s27 = -1, int s28 = -1,
                    int s29 = -1, int s30 = -1, int s31 = -1, int s32 = -1
) {
    return hc_cubic_spline_curve_varg(32, s1, s2, s3, s4,
                s5, s6, s7, s8,
                s9, s10, s11, s12,
                s13, s14, s15, s16,
                s17, s18, s19, s20,
                s21, s22, s23, s24,
                s25, s26, s27, s28,
                s29, s30, s31, s32 );
}

static int hc_lagrange_polynomial_curve(
                    int s1 = -1, int s2 = -1, int s3 = -1, int s4 = -1,
                    int s5 = -1, int s6 = -1, int s7 = -1, int s8 = -1,
                    int s9 = -1, int s10 = -1, int s11 = -1, int s12 = -1,
                    int s13 = -1, int s14 = -1, int s15 = -1, int s16 = -1,
                    int s17 = -1, int s18 = -1, int s19 = -1, int s20 = -1,
                    int s21 = -1, int s22 = -1, int s23 = -1, int s24 = -1,
                    int s25 = -1, int s26 = -1, int s27 = -1, int s28 = -1,
                    int s29 = -1, int s30 = -1, int s31 = -1, int s32 = -1
) {
    return hc_lagrange_curve_varg(32, s1, s2, s3, s4,
                s5, s6, s7, s8,
                s9, s10, s11, s12,
                s13, s14, s15, s16,
                s17, s18, s19, s20,
                s21, s22, s23, s24,
                s25, s26, s27, s28,
                s29, s30, s31, s32 );
}

static int hc_polynomial_curve(
                    double s1 = -1, double s2 = -1, double s3 = -1, double s4 = -1,
                    double s5 = -1, double s6 = -1, double s7 = -1, double s8 = -1,
                    double s9 = -1, double s10 = -1, double s11 = -1, double s12 = -1,
                    double s13 = -1, double s14 = -1, double s15 = -1, double s16 = -1,
                    double s17 = -1, double s18 = -1, double s19 = -1, double s20 = -1,
                    double s21 = -1, double s22 = -1, double s23 = -1, double s24 = -1,
                    double s25 = -1, double s26 = -1, double s27 = -1, double s28 = -1,
                    double s29 = -1, double s30 = -1, double s31 = -1, double s32 = -1
        )

{
        return hc_polynomial_curve_varg(32, s1, s2, s3, s4,
                s5, s6, s7, s8,
                s9, s10, s11, s12,
                s13, s14, s15, s16,
                s17, s18, s19, s20,
                s21, s22, s23, s24,
                s25, s26, s27, s28,
                s29, s30, s31, s32 );

}
static int hc_control_point(double x, double y) {return faust_control_point_map.map(share(control_point{x, y}));}

static int hc_segment(double frac_size, double y_destination, int curve_base_index)
{
    std::cout << "curve base index " << curve_base_index << " & ptr : " << faust_curve_base_map[curve_base_index].get() << std::endl;
    return faust_segment_map.map(share(segment(frac_size, y_destination, faust_curve_base_map[curve_base_index])));
}

// Hypercurve
int hc_curvemaker(int size, double y_start, int n_args, ...)
{
    int i, arg;
    va_list ap;

    std::cout << "number of args : " << n_args << std::endl;
    std::vector<segment > segs;
    va_start(ap, n_args);
    for(i = 0; i < n_args; i++) {
        arg = va_arg(ap, int);
        if(arg <= 0) continue;
        segs.push_back(*faust_segment_map[arg]) ;
    }
    va_end(ap);
    int index = faust_curve_map.map(share(curve(size, y_start, segs)));
    faust_curve_map[index]->ascii_display("faust curve", "curve display", '*');
    return index;
}

static int hc_hypercurve(int size, double y_start,
                    int s1 = -1, int s2 = -1, int s3 = -1, int s4 = -1,
                    int s5 = -1, int s6 = -1, int s7 = -1, int s8 = -1,
                    int s9 = -1, int s10 = -1, int s11 = -1, int s12 = -1,
                    int s13 = -1, int s14 = -1, int s15 = -1, int s16 = -1,
                    int s17 = -1, int s18 = -1, int s19 = -1, int s20 = -1,
                    int s21 = -1, int s22 = -1, int s23 = -1, int s24 = -1,
                    int s25 = -1, int s26 = -1, int s27 = -1, int s28 = -1,
                    int s29 = -1, int s30 = -1, int s31 = -1, int s32 = -1,
                    int s33 = -1, int s34 = -1, int s35 = -1, int s36 = -1,
                    int s37 = -1, int s38 = -1, int s39 = -1, int s40 = -1,
                    int s41 = -1, int s42 = -1, int s43 = -1, int s44 = -1,
                    int s45 = -1, int s46 = -1, int s47 = -1, int s48 = -1,
                    int s49 = -1, int s50 = -1, int s51 = -1, int s52 = -1,
                    int s53 = -1, int s54 = -1, int s55 = -1, int s56 = -1,
                    int s57 = -1, int s58 = -1, int s59 = -1, int s60 = -1,
                    int s61 = -1, int s62 = -1, int s63 = -1, int s64 = -1
) {
    return hc_curvemaker(size, y_start, 64, s1, s2, s3, s4,
                        s5, s6, s7, s8,
                        s9, s10, s11, s12,
                        s13, s14, s15, s16,
                        s17, s18, s19, s20,
                        s21, s22, s23, s24,
                        s25, s26, s27, s28,
                        s29, s30, s31, s32,
                        s33, s34, s35, s36,
                        s37, s38, s39, s40,
                        s41, s42, s43, s44,
                        s45, s46, s47, s48,
                        s49, s50, s51, s52,
                        s53, s54, s55, s56,
                        s57, s58, s59, s60,
                        s61, s62, s63, s64
    );
}

static double hc_run(int curve_index, double read_index)
{
    int i_phasor = std::round(read_index * faust_curve_map[curve_index]->get_definition());
    return faust_curve_map[curve_index]->get_sample_at(i_phasor);
}

static double hc_runi(int index, double phasor)
{
    size_t i_phasor = std::floor(phasor * faust_curve_map[index]->get_definition());
    if( (i_phasor == 0) || (i_phasor >= (faust_curve_map[index]->get_definition() - 1) ))
        return limit(-1, 1, faust_curve_map[index]->get_sample_at(i_phasor));
    size_t n_phasor = i_phasor + 1;
    return limit(-1, 1, linear_interpolation(
             faust_curve_map[index]->get_sample_at(i_phasor),
             faust_curve_map[index]->get_sample_at(n_phasor),
             relative_position(
                 hypercurve::fraction(i_phasor, faust_curve_map[index]->get_definition()),
                 hypercurve::fraction(n_phasor, faust_curve_map[index]->get_definition()),
                 phasor)));

}
}

#endif // HYPERCURVE_FAUST_H
