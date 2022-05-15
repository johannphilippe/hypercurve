/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#include"../src/hypercurve.h"
#include<plugin.h>
#include<OpcodeBase.hpp>
#include<unordered_map>

using namespace hypercurve;

/////////////////////////////////////////////////////////////
// Csound Hypercurve Plugin
/////////////////////////////////////////////////////////////

#include<cstdarg>

// Csound runtime hypercurve
class cs_rt_hypercurve
{
public:

    void _initialize(size_t definition, double y_start_, double *memory)
    {
        _definition = definition;
        y_start = y_start_;
        samples.init(memory, definition);
    }

    void process_init() {it = samples.begin(); }
    void process_one(std::shared_ptr<segment> seg)
    {
        seg->process(it, std::floor(seg->fractional_size * _definition));
    }

    double get_sample_at(size_t i) {return samples[i]; }
    size_t get_definition() {return _definition;}

    void normalize_y(double target_min, double target_max)
    {
        find_extremeness();

        for(size_t i = 0; i < samples.size(); i++)
        {
            samples[i] = ((samples[i] - min ) / ambitus )  * std::abs(target_max - target_min) + target_min;
        }
    }

    void find_extremeness()
    {
        min = samples[0], max = samples[0];
        for(auto & it : samples)
        {
            if(it < min)
                min = it;
            if(it > max)
                max = it;
        }
        ambitus = std::abs(max - min);
    }

protected:

    size_t _definition;
    double y_start;
    memory_vector<double> samples;
    memory_vector<double>::iterator it;
    double min, max, ambitus;
};

// Make it derive from AuxMem
template<typename T>
struct increment_map  : public std::unordered_map<int, T>
{
    increment_map()
    {
    }

    int map( T to_map)
    {
        this->insert({_index, to_map});
        ++_index;
        return _index - 1;
    }
    void unmap(size_t index)
    {
        if(this->find(index) != this->end())
        {
            this->erase(index);
        }
    }

    size_t _index = 0;
};
static increment_map< std::shared_ptr<curve_base> > curve_base_map;
static increment_map< std::shared_ptr<control_point> > control_point_map;
static increment_map< std::shared_ptr<segment> > segment_map;
static increment_map< cs_rt_hypercurve * > curve_map;


//////////////////////////////////////////////////
// Helpers
//////////////////////////////////////////////////
struct cs_control_point : public csnd::Plugin<1, 2>
{
  int init()
  {
      cp = std::make_shared<control_point>(inargs[0], inargs[1]);
      index = control_point_map.map(cp);
      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
      control_point_map.unmap(index);
      return OK;
  }

  int index;
  std::shared_ptr<control_point> cp;
};

//////////////////////////////////////////////////
// Curve lib
//////////////////////////////////////////////////
struct cs_diocles_curve : public csnd::Plugin<1, 1>
{
  int init()
  {
      _curve = std::make_shared<diocles_curve>(inargs[0]);
      index = curve_base_map.map(_curve);
      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
      curve_base_map.unmap(index);
      return OK;
  }

  std::shared_ptr<diocles_curve> _curve;
  int index;
};

struct cs_cubic_curve : public csnd::Plugin<1, 0>
{
    int init()
    {
        _curve = std::make_shared<cubic_curve>();
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<cubic_curve> _curve;
};

struct cs_hanning : public csnd::Plugin<1, 0>
{

    int init()
    {
        _curve = std::make_shared<hanning_curve>();
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<hanning_curve> _curve;
};

struct cs_hamming : public csnd::Plugin<1, 0>
{

    int init()
    {
        _curve = std::make_shared<hamming_curve>();
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<hamming_curve> _curve;
};

struct cs_blackman : public csnd::Plugin<1, 0>
{

    int init()
    {
        _curve = std::make_shared<blackman_curve>();
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<blackman_curve> _curve;
};

struct cs_gaussian : csnd::Plugin<1,2>
{
    int init()
    {
        _curve = std::make_shared<gauss_curve>(inargs[0], inargs[1]);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<gauss_curve> _curve;
};

struct cs_typed : csnd::Plugin<1,1>
{
    int init()
    {
        _curve = std::make_shared<typed_curve>(inargs[0]);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<typed_curve> _curve;
};

struct cs_power : csnd::Plugin<1,1>
{
    int init()
    {
        _curve = std::make_shared<power_curve>(inargs[0]);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<power_curve> _curve;
};

struct cs_toxoid : csnd::Plugin<1, 1>
{
    int init()
    {
        _curve = std::make_shared<toxoid_curve>(inargs[0]);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<toxoid_curve> _curve;
};


struct cs_mouse : csnd::Plugin<1, 0>
{
    int init()
    {
        _curve = std::make_shared<mouse_curve>();
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<mouse_curve> _curve;
};

struct cs_bicorn : csnd::Plugin<1, 1>
{
    int init()
    {
        _curve = std::make_shared<bicorn_curve>(inargs[0] == 1);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<bicorn_curve> _curve;
};

struct cs_catenary : csnd::Plugin<1, 1>
{
    int init()
    {
        _curve = std::make_shared<catenary_curve>(inargs[0]);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<catenary_curve> _curve;
};

struct cs_tightrope_walker : csnd::Plugin<1, 2>
{
    int init()
    {
        _curve = std::make_shared<tightrope_walker_curve>(inargs[0], inargs[1]);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<tightrope_walker_curve> _curve;
};

struct cs_quad_bezier : csnd::Plugin<1, 1>
{
    int init()
    {
        cp1 = control_point_map[inargs[0]];
        _curve = std::make_shared<quadratic_bezier_curve>(*cp1);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<control_point> cp1;
    std::shared_ptr<quadratic_bezier_curve> _curve;
};

struct cs_cub_bezier : csnd::Plugin<1,2>
{
    int init()
    {
        cp1 = control_point_map[inargs[0]];
        cp2 = control_point_map[inargs[1]];
        _curve = std::make_shared<cubic_bezier_curve>(*cp1, *cp2);
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<control_point> cp1, cp2;
    std::shared_ptr<cubic_bezier_curve> _curve;
};

struct cs_catmull_rom : csnd::Plugin<1, 2>
{

    int init()
    {
        cp1 = control_point_map[inargs[0]];
        cp2 = control_point_map[inargs[1]];
        _curve = std::make_shared<catmull_rom_spline_curve>(*cp1, *cp2);
        index= curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    std::shared_ptr<control_point> cp1, cp2;
    std::shared_ptr<catmull_rom_spline_curve> _curve;
};

struct cs_polynomial_curve : csnd::Plugin<1, 64>
{
    int init()
    {
        args.allocate(csound, in_count());
        std::cout << "in count polynomial : " << in_count() << std::endl;
        for(size_t i = 0; i < args.len(); ++i)
            args[i] = inargs[i];
        _curve = std::make_shared<polynomial_curve>(args.data(), args.len());
        index = curve_base_map.map(_curve);
        outargs[0] = index;
        return OK;
    }

    int deinit()
    {
        curve_base_map.unmap(index);
        return OK;
    }

    int index;
    csnd::AuxMem<double> args;
    std::shared_ptr<polynomial_curve> _curve;
};

//////////////////////////////////////////////////
// Invert curve
//////////////////////////////////////////////////

struct cs_invert_curve : public csnd::Plugin<1,1>
{
  int init()
  {
      if(curve_base_map.find(int(inargs[0])) == curve_base_map.end())
          return NOTOK;
      std::shared_ptr<curve_base> c = curve_base_map[int(inargs[0])];
      c->inverted = true;
      outargs[0] = int(inargs[0]);
      return OK;
  }
};

//////////////////////////////////////////////////
// Normalize curve
//////////////////////////////////////////////////

struct cs_normalize_curve : public csnd::InPlug< 3>
{
  int init()
  {
      if(curve_map.find(int(args[0])) == curve_map.end())
      {
          return NOTOK;
      }
      curve_map[int(args[0])]->normalize_y(args[1], args[2]);
      return OK;
  }
};



//////////////////////////////////////////////////
// Curve lib
//////////////////////////////////////////////////

struct cs_segment : public csnd::Plugin<1, 3>
{
  int init()
  {
      seg = std::make_shared<segment>(inargs[0], inargs[1], curve_base_map[inargs[2]]);
      index = segment_map.map(seg);
      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
      segment_map.unmap(index);
      return OK;
  }

private:
  std::shared_ptr<segment> seg;
  int index;
};

struct cs_curve : public csnd::Plugin<1,64>, public cs_rt_hypercurve
{
  int init()
  {
      _mem.allocate(csound, (int)inargs[0]);
      _initialize(inargs[0], inargs[1], _mem.data());
      check_total_size();
      process_init();


      for(size_t i = 2; i < in_count(); ++i)
      {
          if(segment_map.find(int(inargs[i])) == segment_map.end())
              return NOTOK;

          y_start = (double) ((i > 2) ?  segment_map[inargs[i - 1] ]->y_destination : inargs[1]);
          segment_map[inargs[i]]->init(y_start, segment_map[inargs[i]]->fractional_size * _definition);
          process_one(segment_map[inargs[i]]);
      }
      index = curve_map.map(this);

      // Weirdly, it works below
      AsciiPlotter p("hc_hypercurve : "  + std::to_string(index)  , 80, 15);
      p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()), "helloplot", '*');
      p.show();

      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
      curve_map.unmap(index);
      return OK;
  }

  int index;
  double y_start;
  csnd::AuxMem<double> _mem;

protected:

    void check_total_size()
    {
        double x = 0;
        for(size_t i = 2; i < in_count(); ++i)
            x += segment_map[inargs[i]]->fractional_size;
        if( (x > 1.0) ||  (x < 1.0) )
        {
            this->rescale(x);
        }
    }

    void rescale(double x)
    {
        double factor = 1. / x;
        for(size_t i = 2; i < in_count(); ++i) {
            segment_map[inargs[i]]->rescale_x(factor);
        }
    }
};

// Do the same as below with trigger (trig_run_hypercurve) with frequency and trigger parameter.
struct cs_run_curve : public csnd::Plugin<1, 2>
{
  int init()
  {
      index = inargs[0];
      if(curve_map.find(index) == curve_map.end())
          return NOTOK;
      return OK;
  }

  int kperf()
  {
      phasor = inargs[1];
      i_phasor = std::round(phasor * curve_map[index]->get_definition());
      outargs[0] = limit(-1, 1, curve_map[index]->get_sample_at(i_phasor));
      return OK;
  }

  int aperf()
  {
      for(size_t i = 0; i < nsmps; ++i)
      {
          phasor = inargs(1)[i];
          i_phasor = std::round(phasor * curve_map[index]->get_definition());
          outargs(0)[i] = limit(-1, 1, curve_map[index]->get_sample_at(i_phasor));

      }

      return OK;
  }

private:
  double phasor;
  size_t i_phasor;
  int index;
};


struct cs_run_curve_interpolate : public csnd::Plugin<1, 2>
{
  int init()
  {
      index = inargs[0];
      if(curve_map.find(index) == curve_map.end())
          return NOTOK;
      return OK;
  }

  int kperf()
  {
      phasor = inargs[1];
      outargs[0] = interpolate();
      return OK;
  }

  int aperf()
  {
      for(size_t i = 0; i < nsmps; ++i)
      {
          phasor = inargs(1)[i];
          outargs(0)[i] = interpolate();
      }

      return OK;
  }

  double interpolate()
  {
      i_phasor = floor(phasor * curve_map[index]->get_definition());
      if( (phasor == 0) || (i_phasor >= (curve_map[index]->get_definition() - 1) ))
          return limit(-1, 1, curve_map[index]->get_sample_at(i_phasor));
      n_phasor = ceil(phasor * curve_map[index]->get_definition());
      return limit(-1, 1, linear_interpolation(
                       curve_map[index]->get_sample_at(i_phasor),
                       curve_map[index]->get_sample_at(n_phasor),
                       relative_position(
                           fraction(i_phasor, curve_map[index]->get_definition()),
                           fraction(n_phasor, curve_map[index]->get_definition()),
                           phasor)));
  }

private:
  double phasor;
  size_t i_phasor, n_phasor;
  int index;
};

////////////////////////////////////////////////////////////////////////
// Operators for curves
////////////////////////////////////////////////////////////////////////

struct cs_operator : public cs_rt_hypercurve
{
    void _allocate( csnd::Csound *cs)
    {
        mem.allocate(cs, _definition);
        samples.init(mem.data(), _definition);
    }

    csnd::AuxMem<double> mem;
};

struct cs_add_curve : public csnd::Plugin<1, 2>, cs_operator
{
    int init()
    {
      if( (curve_map.find(int(inargs[0])) == curve_map.end())
              || (curve_map.find(int(inargs[1])) == curve_map.end())
              )
          return NOTOK;
      cs_rt_hypercurve *c1 = curve_map[int(inargs[0])];
      cs_rt_hypercurve *c2 = curve_map[int(inargs[1])];

      _definition = c1->get_definition();
      _allocate(csound);

      for(size_t i = 0; i < _definition; ++i) {
          samples[i] = c1->get_sample_at(i) + c2->get_sample_at(i);
      }
      index = curve_map.map(this);

      AsciiPlotter p("hc_hypercurve : "  + std::to_string(index)  , 80, 15);
      p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()), "helloplot", '*');
      p.show();

      outargs[0] = index;
      return OK;
    }

    int index;
};

struct cs_sub_curve : public csnd::Plugin<1, 2>, cs_operator
{
    int init()
    {
      if( (curve_map.find(int(inargs[0])) == curve_map.end())
              || (curve_map.find(int(inargs[1])) == curve_map.end())
              )
          return NOTOK;
      cs_rt_hypercurve *c1 = curve_map[int(inargs[0])];
      cs_rt_hypercurve *c2 = curve_map[int(inargs[1])];

      _definition = c1->get_definition();
      _allocate(csound);

      for(size_t i = 0; i < _definition; ++i) {
          samples[i] = c1->get_sample_at(i) - c2->get_sample_at(i);
      }
      index = curve_map.map(this);

      AsciiPlotter p("hc_hypercurve : "  + std::to_string(index)  , 80, 15);
      p.addPlot(std::vector<double>(this->samples.begin(), this->samples.end()), "helloplot", '*');
      p.show();

      outargs[0] = index;
      return OK;
    }

    int index;
    csnd::AuxMem<double> mem;
};

struct cs_mult_curve : public csnd::Plugin<1, 2>, cs_operator
{
    int init()
    {
      if( (curve_map.find(int(inargs[0])) == curve_map.end())
              || (curve_map.find(int(inargs[1])) == curve_map.end())
              )
          return NOTOK;
      cs_rt_hypercurve *c1 = curve_map[int(inargs[0])];
      cs_rt_hypercurve *c2 = curve_map[int(inargs[1])];

      _definition = c1->get_definition();
      _allocate(csound);

      for(size_t i = 0; i < _definition; ++i) {
          samples[i] = c1->get_sample_at(i) * c2->get_sample_at(i);
      }
      index = curve_map.map(this);

      AsciiPlotter p("hc_hypercurve : "  + std::to_string(index)  , 80, 15);
      p.addPlot(std::vector<double>(this->samples.begin(), this->samples.end()), "helloplot", '*');
      p.show();

      outargs[0] = index;
      return OK;
    }

    int index;
    csnd::AuxMem<double> mem;
};

struct cs_div_curve : public csnd::Plugin<1, 2>, cs_operator
{
    int init()
    {
      if( (curve_map.find(int(inargs[0])) == curve_map.end())
              || (curve_map.find(int(inargs[1])) == curve_map.end())
              )
          return NOTOK;
      cs_rt_hypercurve *c1 = curve_map[int(inargs[0])];
      cs_rt_hypercurve *c2 = curve_map[int(inargs[1])];

      _definition = c1->get_definition();
      _allocate(csound);

      for(size_t i = 0; i < _definition; ++i) {
          samples[i] = c1->get_sample_at(i) / c2->get_sample_at(i);
      }
      index = curve_map.map(this);

      AsciiPlotter p("hc_hypercurve : "  + std::to_string(index)  , 80, 15);
      p.addPlot(std::vector<double>(this->samples.begin(), this->samples.end()), "helloplot", '*');
      p.show();

      outargs[0] = index;
      return OK;
    }

    int index;
    csnd::AuxMem<double> mem;
};

#include <modload.h>

void csnd::on_load(Csound *csound) {
    std::cout << "loading csound hypercurve" << std::endl;
    // Helpers
    csnd::plugin<cs_normalize_curve>(csound, "hc_normalize_y", "", "iii", csnd::thread::i);
    csnd::plugin<cs_control_point>(csound, "hc_control_point", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_control_point>(csound, "hc_point", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_invert_curve>(csound, "hc_invert", "i", "i", csnd::thread::i);
    csnd::plugin<cs_add_curve>(csound, "hc_add", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_sub_curve>(csound, "hc_sub", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_mult_curve>(csound, "hc_mult", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_div_curve>(csound, "hc_div", "i", "ii" , csnd::thread::i);

    // Curve types
    csnd::plugin<cs_diocles_curve>(csound, "hc_diocles_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_cubic_curve>(csound, "hc_cubic_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_power>(csound, "hc_power_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_hanning>(csound, "hc_hanning_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_hamming>(csound, "hc_hamming_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_blackman>(csound, "hc_blackman_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_gaussian>(csound, "hc_gauss_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_catenary>(csound, "hc_catenary_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_toxoid>(csound, "hc_toxoid_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_mouse>(csound, "hc_mouse_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_bicorn>(csound, "hc_bicorn_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_tightrope_walker>(csound, "hc_tightrope_walker_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_polynomial_curve>(csound, "hc_polynomial_curve", "i", "m", csnd::thread::i);

    csnd::plugin<cs_typed>(csound, "hc_typed_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_quad_bezier>(csound, "hc_quadratic_bezier_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_cub_bezier>(csound, "hc_cubic_bezier_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_catmull_rom>(csound, "hc_catmull_rom_curve", "i" , "ii", csnd::thread::i);

    // Aliases
    csnd::plugin<cs_toxoid>(csound, "hc_duplicatrix_cubic_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_catenary>(csound, "hc_funicular_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_diocles_curve>(csound, "hc_cissoid_curve", "i", "i", csnd::thread::i); // alias for diocles
    csnd::plugin<cs_gaussian>(csound, "hc_gaussian_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_typed>(csound, "hc_transeg_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_mouse>(csound, "hc_kiss_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_bicorn>(csound, "hc_cocked_hat_curve", "i", "i", csnd::thread::i);

    // Core
    csnd::plugin<cs_segment>(csound, "hc_segment", "i", "iii", csnd::thread::i);
    csnd::plugin<cs_curve>(csound, "hc_hypercurve", "i", "iim", csnd::thread::i);

    // Run
    csnd::plugin<cs_run_curve>(csound, "hc_run_hypercurve", "k", "ik", csnd::thread::ik);
    csnd::plugin<cs_run_curve>(csound, "hc_run", "k", "ik", csnd::thread::ik);
    csnd::plugin<cs_run_curve>(csound, "hc_run_hypercurve", "a", "ia", csnd::thread::ia);
    csnd::plugin<cs_run_curve>(csound, "hc_run", "a", "ia", csnd::thread::ia);
    // With interpolation
    csnd::plugin<cs_run_curve_interpolate>(csound, "hc_runi", "k", "ik", csnd::thread::ik);
    csnd::plugin<cs_run_curve_interpolate>(csound, "hc_runi", "a", "ia", csnd::thread::ia);
}
