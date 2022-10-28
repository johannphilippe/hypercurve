/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#include"../src/hypercurve.h"
#include<plugin.h>
#include<OpcodeBase.hpp>

using namespace hypercurve;

/////////////////////////////////////////////////////////////
// Csound Hypercurve Plugin
/////////////////////////////////////////////////////////////

#include<cstdarg>


// Custom FTable utilities
struct ftable_info
{
    int res; // Result for GetTable
    int fno;
    double *t_ptr;
    size_t size = 0;
};

ftable_info allocate_gen(CSOUND_ *csound, int fno, size_t size)
{
    int t_nbr;
    double t_dbl;
    double *t_ptr = nullptr;
    if(fno > 0) {
        t_nbr = fno;
    } else {
        // Automatic table number
        t_nbr = 100;
        do {
            ++t_nbr;
            t_dbl = t_nbr;
            csound->GetTable(csound, &t_ptr, t_nbr);
            FUNC* t = csound->FTnp2Finde(csound, &t_dbl);
            if(t == NULL)
                t_ptr = NULL;

        } while(t_ptr != NULL);
    }


    int res = csound->FTAlloc(csound,  t_nbr, size);
    if(res != 0) return {csound->InitError(csound, "Error alloc Hypercurve GEN %d with res = %d", t_nbr, res ) , fno, nullptr, 0};
    size_t sz = csound->GetTable(csound, &t_ptr, t_nbr );
    return {0, t_nbr, t_ptr, sz};
}

int free_gen(CSOUND_ *csound, int fno)
{
    return csound->FTDelete(csound, fno);
}

ftable_info get_gen(CSOUND_ *csound, int fno)
{
    double *t_ptr = nullptr;
    size_t size = csound->GetTable(csound, &t_ptr, fno);
    if(t_ptr == (MYFLT*)NULL)
    {
        return {
            csound->InitError(csound, "Error init Hypercurve GEN %d", fno ) , fno, nullptr, 0
        };
    }
    return {0, fno, t_ptr, size};
}


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
        //
        seg->process(it, std::round(seg->fractional_size * _definition));
    }

    double *get_samples() {return samples.data();}
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

    ftable_info table;
protected:

    size_t _definition;
    double y_start;
    memory_vector<double> samples;
    memory_vector<double>::iterator it;
    double min, max, ambitus;
};

static increment_map< std::shared_ptr<curve_base> > curve_base_map;
static increment_map< std::shared_ptr<control_point> > control_point_map;
static increment_map< std::shared_ptr<segment> > segment_map;
static std::unordered_map< int,  cs_rt_hypercurve * > curve_map;

//////////////////////////////////////////////////
// Helpers
//////////////////////////////////////////////////
struct cs_control_point : public csnd::Plugin<1, 2>
{
  int init()
  {
      csound->plugin_deinit(this);
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

struct cs_linear_curve : public csnd::Plugin<1, 0>
{
  int init()
  {
      csound->plugin_deinit(this);
      _curve = std::make_shared<linear_curve>();
      index = curve_base_map.map(_curve);
      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
    curve_base_map.unmap(index);
    return OK;
  }

  std::shared_ptr<linear_curve> _curve;
  int index;
};

struct cs_diocles_curve : public csnd::Plugin<1, 1>
{
  int init()
  {
      csound->plugin_deinit(this);
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


struct cs_logarithmic_curve : public csnd::Plugin<1, 0>
{
  int init()
  {
      csound->plugin_deinit(this);
      _curve = std::make_shared<logarithmic_curve>();
      index = curve_base_map.map(_curve);
      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
      curve_base_map.unmap(index);
      return OK;
  }

  std::shared_ptr<logarithmic_curve> _curve;
  int index;
};

struct cs_exponential_curve : public csnd::Plugin<1, 0>
{
  int init()
  {
      csound->plugin_deinit(this);
      _curve = std::make_shared<exponential_curve>();
      index = curve_base_map.map(_curve);
      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
      curve_base_map.unmap(index);
      return OK;
  }

  std::shared_ptr<exponential_curve> _curve;
  int index;
};

struct cs_cubic_curve : public csnd::Plugin<1, 0>
{
    int init()
    {
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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


struct cs_mouth : csnd::Plugin<1, 0>
{
    int init()
    {
        csound->plugin_deinit(this);
        _curve = std::make_shared<mouth_curve>();
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
    std::shared_ptr<mouth_curve> _curve;
};

struct cs_bicorn : csnd::Plugin<1, 1>
{
    int init()
    {
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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
        csound->plugin_deinit(this);
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

struct cs_lagrange_curve : csnd::Plugin<1, 64>
{
    int init()
    {
        csound->plugin_deinit(this);
        cp.allocate(csound, in_count() + 2);
        for(size_t i = 0; i < in_count(); ++i)
            cp[i+1] = *control_point_map[inargs[i]];
        _curve = std::make_shared<lagrange_polynomial_curve>(cp.data(), cp.len());
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
    csnd::AuxMem<control_point> cp;
    std::shared_ptr<lagrange_polynomial_curve> _curve;
};

struct cs_cubic_spline_curve : csnd::Plugin<1, 64>
{
    int init()
    {
        csound->plugin_deinit(this);
        int n = in_count() + 2 - 1;
        int nn = (5 * n) + (3 * ( n + 1 )) + (n - 1)  + ( (n - 1) * n);
        spline_memory.allocate(csound, nn);

        cp.allocate(csound, in_count() + 2);
        for(size_t i = 0; i < in_count(); ++i)
            cp[i+1] = *control_point_map[inargs[i]];

        _curve = std::make_shared<cubic_spline_curve>(cp.data(), cp.len(), spline_memory.data(), spline_memory.len());
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
    csnd::AuxMem<double> spline_memory;
    csnd::AuxMem<control_point> cp;
    std::shared_ptr<cubic_spline_curve> _curve;
};

struct cs_polynomial_curve : csnd::Plugin<1, 64>
{
    int init()
    {
        csound->plugin_deinit(this);
        args.allocate(csound, in_count());
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

struct cs_mirror_curve : public csnd::Plugin<1, 1>
{
    int init()
    {
      if(curve_base_map.find(int(inargs[0])) == curve_base_map.end())
          return NOTOK;
      std::shared_ptr<curve_base> c = curve_base_map[int(inargs[0])];
      c->mirrored = !c->mirrored;
      outargs[0] = int(inargs[0]);
      return OK;
    }
};

//////////////////////////////////////////////////
// Normalize curve
//////////////////////////////////////////////////

struct cs_normalize_curve : public csnd::InPlug<3>
{
  int init()
  {
      if(curve_map.find(int(args[0])) == curve_map.end())
      {
          std::cout << "Normalize cannot find the expected curve " << std::endl;
          return NOTOK;
      }
      cs_rt_hypercurve *crv = curve_map[int(args[0])];
      crv->normalize_y(args[1], args[2]);
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
      csound->plugin_deinit(this);
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

// Hypercurve implemented as GEN
struct cs_gen : public csnd::Plugin<1, 67>, public cs_rt_hypercurve
{
  int init()
  {
      csound->plugin_deinit(this);
      table = allocate_gen(csound->get_csound(), inargs[0], inargs[1]);
      std::cout << "gen allocated with result : " << table.res << std::endl;
      if(table.res != 0) {
          return table.res;
      }


      _initialize(inargs[1], inargs[2], table.t_ptr);
      check_total_size();
      process_init();

      for(size_t i = 3; i < in_count(); ++i)
      {
          if(segment_map.find(int(inargs[i])) == segment_map.end())
              return NOTOK;

          y_start = (double) ((i > 3) ?  segment_map[inargs[i - 1] ]->y_destination : inargs[2]);
          segment_map[inargs[i]]->init(y_start, segment_map[inargs[i]]->fractional_size * _definition);
          process_one(segment_map[inargs[i]]);
      }

      // Weirdly, it works below
      const char * name = csound->get_csound()->GetOutputArgName(this, 0);
      AsciiPlotter p(std::string(name)  +  " - Hypercurve GEN number : "  + std::to_string(table.fno)  , 80, 15);
      p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()), std::string(" - ") + name , (char)(rand() % (126 - 92) + 92));
      p.show();

      curve_map[table.fno] = this;
      outargs[0] = table.fno;
      return OK;
  }

    void check_total_size()
    {
        double x = 0;
        for(size_t i = 3; i < in_count(); ++i)
            x += segment_map[inargs[i]]->fractional_size;
        if( (x > 1.0) ||  (x < 1.0) )
        {
            this->rescale(x);
        }
    }

    void rescale(double x)
    {
        double factor = 1. / x;
        for(size_t i = 3; i < in_count(); ++i) {
            segment_map[inargs[i]]->rescale_x(factor);
        }
    }

    int deinit()
    {
        std::cout << "hypercurve deallocation for gen number : " << table.fno << std::endl;
        csound->get_csound()->FTDelete(csound->get_csound(), table.fno);
        curve_map.erase(table.fno);
        return OK;
    }

};


struct cs_concat : public csnd::Plugin<1, 66>, public cs_rt_hypercurve
{
    int init()
    {
        csound->plugin_deinit(this);
        size_t total_size = count_total_size();
        size_t size = (inargs[1] == 0) ? total_size : inargs[1];
        table = allocate_gen(csound->get_csound(), inargs[0], size);
        std::cout << "gen allocated with result : " << table.res << std::endl;
        if(table.res != 0) {
            return table.res;
        }

        t = get_gen(csound->get_csound(), inargs[2]);
        _initialize(size, *t.t_ptr, table.t_ptr);

        if(inargs[1] == 0)
            concat_noresize();
        else
            concat_resize(total_size);

        // Weirdly, it works below
        const char * name = csound->get_csound()->GetOutputArgName(this, 0);
        AsciiPlotter p(std::string(name)  +  " - Hypercurve GEN number : "  + std::to_string(table.fno)  , 80, 15);
        p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()), std::string(" - ") + name , (char)(rand() % (126 - 92) + 92));
        p.show();

        curve_map[table.fno] = this;
        outargs[0] = table.fno;
        return OK;
    }


    void concat_resize(size_t total_size)
    {
        size_t current = 2;
        double incr = total_size / inargs[1];
        int passed = 0;
        for(size_t i = 0; i < table.size; ++i)
        {

            double index = i * incr;
            double cindex = index - passed;
            t = get_gen(csound->get_csound(), inargs[current]);
            if(cindex >= t.size)
            {
                current++;
                passed += t.size;
                t = get_gen(csound->get_csound(), inargs[current]);
            }

            double value;
            if(cindex == 0)
            {
                value = *t.t_ptr;
            } else {
                int f = floor(cindex);
                int c = f + 1;
                double relative_x = relative_position(f, c, cindex);
                value = linear_interpolation(t.t_ptr[f], t.t_ptr[c], relative_x);
            }
            this->samples[i] = value;
        }

    }

    void concat_noresize()
    {
        memory_vector<double>::iterator it = samples.begin();
        for(size_t i = 2; i < in_count(); ++i)
        {
                t = get_gen(csound->get_csound(), inargs[i]);
            for(size_t j = 0; j < t.size; ++j)
            {
                *it = t.t_ptr[j];
                ++it;
            }
        }
    }

    size_t count_total_size()
    {
        size_t sum = 0;
        for(size_t i = 2; i < in_count(); ++i)
        {
            ftable_info t = get_gen(csound->get_csound(), inargs[i]);
            sum += t.size;
        }
        return sum;
    }

    int deinit()
    {
        std::cout << "hypercurve deallocation for gen number : " << table.fno << std::endl;
        csound->get_csound()->FTDelete(csound->get_csound(), table.fno);
        curve_map.erase(table.fno);
        return OK;
    }

    ftable_info t;
};

// Modify with temp buffer so that it does not create a new hypercurve
struct cs_resize : public csnd::Plugin<1, 3>, public cs_rt_hypercurve
{

    int init()
    {
        table = allocate_gen(csound->get_csound(), inargs[0], inargs[1]);
        tt = get_gen(csound->get_csound(), inargs[2]);
        _initialize(inargs[1], *tt.t_ptr, this->table.t_ptr);
        double incr = tt.size / inargs[1];
        for(size_t i = 0; i < inargs[1]; ++i)
        {
            double index = double(i) * incr;

            int f = floor(index), c = ceil(index);
            double relative_x = relative_position(f, c, index);
            double value = linear_interpolation(tt.t_ptr[f], tt.t_ptr[c], relative_x );
            this->samples[i] = value;
        }
        // Weirdly, it works below
        const char * name = csound->get_csound()->GetOutputArgName(this, 0);
        AsciiPlotter p(std::string(name)  +  " - Hypercurve GEN number : "  + std::to_string(table.fno)  , 80, 15);
        p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()), std::string(" - ") + name , (char)(rand() % (126 - 92) + 92));
        p.show();

        outargs[0] = table.fno;
        return OK;
    }

    int deinit()
    {
        csound->get_csound()->FTDelete(csound->get_csound(), table.fno);
        return OK;
    }
    ftable_info tt;
};


////////////////////////////////////////////////////////////////////////
// Operators for curves
////////////////////////////////////////////////////////////////////////

struct cs_operator_add { double process(double a, double b) {return a + b;}};
struct cs_operator_sub { double process(double a, double b) {return a - b;}};
struct cs_operator_mult { double process(double a, double b) {return a * b;}};
struct cs_operator_div { double process(double a, double b) {return a / b;}};

template<typename T = cs_operator_add>
struct cs_operator
        : public csnd::Plugin<1, 2>
        , public cs_rt_hypercurve
        , public T
{
    int init()
    {
        t1 = get_gen(csound->get_csound(), inargs[0]);
        t2 = get_gen(csound->get_csound(), inargs[1]);
      if(t1.res != 0 || t2.res != 0)
          return NOTOK;

      _definition = std::max(t1.size, t2.size);

      table = allocate_gen(csound->get_csound(), 0, _definition);
      if(table.res != 0)
          return NOTOK;

      _initialize(_definition, y_start, table.t_ptr);
      for(size_t i = 0; i < table.size; ++i) {
          table.t_ptr[i] = this->process(t1.t_ptr[i], t2.t_ptr[i]);
      }

      // Weirdly, it works below
      const char * name = csound->get_csound()->GetOutputArgName(this, 0);
      AsciiPlotter p(std::string(name)
                     +  " - Hypercurve GEN number : "
                     +  std::to_string(table.fno)  , 80, 15);
      p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()),
                 std::string(" - ") + name ,
                 (char)(rand() % (126 - 92) + 92));
      p.show();

      curve_map[table.fno] = this;
      outargs[0] = table.fno;
      return OK;

    }
    ftable_info t1, t2;
};

template<typename T = cs_operator_add>
struct cs_operator_num
        : public csnd::Plugin<1, 2>
        , public cs_rt_hypercurve
        , public T
{
    int init()
    {
        t1 = get_gen(csound->get_csound(), inargs[0]);
        double t2 = inargs[1];
      if(t1.res != 0 )
          return NOTOK;

      _definition = t1.size;

      table = allocate_gen(csound->get_csound(), 0, _definition);
      if(table.res != 0)
          return NOTOK;

      _initialize(_definition, y_start, table.t_ptr);
      for(size_t i = 0; i < table.size; ++i) {
          table.t_ptr[i] = this->process(t1.t_ptr[i], t2);
      }

      // Weirdly, it works below
      const char * name = csound->get_csound()->GetOutputArgName(this, 0);
      AsciiPlotter p(std::string(name)
                     +  " - Hypercurve GEN number : "
                     +  std::to_string(table.fno)  , 80, 15);
      p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()),
                 std::string(" - ") + name ,
                 (char)(rand() % (126 - 92) + 92));
      p.show();

      curve_map[table.fno] = this;
      outargs[0] = table.fno;
      return OK;

    }
    ftable_info t1, t2;
};

struct cs_add_curve : public cs_operator<cs_operator_add>{};
struct cs_sub_curve : public cs_operator<cs_operator_sub>{};
struct cs_mult_curve : public cs_operator<cs_operator_mult>{};
struct cs_div_curve : public cs_operator<cs_operator_div>{};
struct cs_add_curve_num : public cs_operator_num<cs_operator_add>{};
struct cs_sub_curve_num : public cs_operator_num<cs_operator_sub>{};
struct cs_mult_curve_num : public cs_operator_num<cs_operator_mult>{};
struct cs_div_curve_num : public cs_operator_num<cs_operator_div>{};

// args : crv, path, waveform , fill, grid, invert
struct cs_export_png : public csnd::InPlug<6>
{
    int init()
    {
        if(curve_map.find(int(args[0])) == curve_map.end()) {
            return NOTOK;
        }
        cs_rt_hypercurve *crv = curve_map[int(args[0])];

        std::string path(args.str_data(1).data);
        hypercurve::png png(2048, 1024, args[5] != 0 ? white : black, args[5] != 0 ? red : purple);

        png.draw_curve(crv->table.t_ptr, crv->table.size, int(args[3]) != 0, int(args[2]) != 0);
        if(args[4] != 0) {
            png.draw_grid(10, 10, args[5] != 0 ? black : white);
        }
        png.write_as_png(path);
        return OK;
    }
};

// cs_spat for several channel
/*
 * arg1 = number of channels
 * arg2 = function table
 * arg3 = ktrajectory (0-1 range)
*/

/*
struct cs_spat : public csnd::Plugin<1, 3>
{
  int init()
  {

      return OK;
  }
};

*/

// fno, size, max_segs, min, max, is_envelope, is_waveform
struct cs_random_curve_composer : public csnd::Plugin<1, 7>,  public cs_rt_hypercurve
{
    int init()
    {
        table = allocate_gen(csound->get_csound(), inargs[0], inargs[1]);
        if(table.res != 0) {
            return NOTOK;
        }

        _initialize(inargs[1], inargs[2], table.t_ptr);
        process_init();

        int max_segs = inargs[2];
        const double min = inargs[3];
        const double max = inargs[4];
        const bool envelop = inargs[5] == 1;
        const bool waveform = inargs[6] == 1;

        std::string cnames = "";
        size_t nsegs = 1 + (rand() % max_segs);
        if(nsegs < 2) nsegs = 2;

        y_start = (envelop || waveform) ? 0 : random<double>(0, 1);
        for(size_t i = 0; i < nsegs; ++i) {
            double frac_size = random<double>(0.1, 1);
            double dest = waveform ? random<double>(-1, 1) : random<double>(0, 1);
            curve_base_index index =  gen_curve();

            while( index == user_defined_i ||
                   (index == polynomial_i) ||
                   (index == lagrange_polynomial_i) ||
                   (index == cubic_spline_i) ||
                   index == catmull_rom_spline_i)
            {
                index = gen_curve();
            }


            cnames += std::string(get_curve_base_name(index)) + "_";
            std::pair<std::vector<double>, std::vector<control_point>> args = random_args_generator(index);
            if(i == (nsegs - 1) && (envelop | waveform) )
                dest = 0;

            auto seg = std::make_shared<segment>(frac_size, dest, get_curve_from_index(index, args.first, args.second) );
            process_one(seg);
        }
        normalize_y(min, max);

        const char * name = csound->get_csound()->GetOutputArgName(this, 0);
        AsciiPlotter p(std::string(name)  +  " - Hypercurve GEN number : "  + std::to_string(table.fno)  , 80, 15);
        p.addPlot( std::vector<double>(this->samples.begin(), this->samples.end()),
                   std::string(" - ") + name ,
                   (char)(rand() % (126 - 92) + 92));
        p.show();
        outargs[0] = table.fno;
        return OK;
    }


    inline curve_base_index gen_curve()
    {
        return static_cast<curve_base_index>(rand() % static_cast<int>(size_i));
    }

};

#include <modload.h>

void csnd::on_load(Csound *csound) {
    std::cout << "💫 Loading HYPERCURVE 💫" << std::endl;
    // Helpers
    csnd::plugin<cs_normalize_curve>(csound, "hc_normalize", "", "iii", csnd::thread::i);

    csnd::plugin<cs_control_point>(csound, "hc_control_point", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_control_point>(csound, "hc_point", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_invert_curve>(csound, "hc_invert", "i", "i", csnd::thread::i);
    csnd::plugin<cs_mirror_curve>(csound, "hc_mirror", "i", "i", csnd::thread::i);
    csnd::plugin<cs_add_curve>(csound, "hc_add", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_sub_curve>(csound, "hc_sub", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_mult_curve>(csound, "hc_mult", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_div_curve>(csound, "hc_div", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_add_curve_num>(csound, "hc_addn", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_sub_curve_num>(csound, "hc_subn", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_mult_curve_num>(csound, "hc_multn", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_div_curve_num>(csound, "hc_divn", "i", "ii" , csnd::thread::i);
    csnd::plugin<cs_concat>(csound, "hc_concat", "i", "iim" , csnd::thread::i);


    csnd::plugin<cs_export_png>(csound, "hc_write_as_png", "", "iSoppo", csnd::thread::i);
    //csnd::plugin<cs_random_curve_composer>(csound, "hc_random_curve", "i", "iiiiioo", csnd::thread::i);

    // Curve types
    csnd::plugin<cs_diocles_curve>(csound, "hc_diocles_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_logarithmic_curve>(csound, "hc_logarithmic_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_exponential_curve>(csound, "hc_exponential_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_linear_curve>(csound, "hc_linear_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_cubic_curve>(csound, "hc_cubic_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_power>(csound, "hc_power_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_hanning>(csound, "hc_hanning_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_hamming>(csound, "hc_hamming_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_blackman>(csound, "hc_blackman_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_gaussian>(csound, "hc_gauss_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_catenary>(csound, "hc_catenary_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_toxoid>(csound, "hc_toxoid_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_mouth>(csound, "hc_mouth_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_bicorn>(csound, "hc_bicorn_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_tightrope_walker>(csound, "hc_tightrope_walker_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_polynomial_curve>(csound, "hc_polynomial_curve", "i", "m", csnd::thread::i);

    csnd::plugin<cs_typed>(csound, "hc_typed_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_quad_bezier>(csound, "hc_quadratic_bezier_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_cub_bezier>(csound, "hc_cubic_bezier_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_catmull_rom>(csound, "hc_catmull_rom_curve", "i" , "ii", csnd::thread::i);
    csnd::plugin<cs_lagrange_curve>(csound, "hc_lagrange_curve", "i", "m", csnd::thread::i);
    csnd::plugin<cs_cubic_spline_curve>(csound, "hc_cubic_spline_curve", "i", "m", csnd::thread::i);

    // Aliases
    csnd::plugin<cs_toxoid>(csound, "hc_duplicatrix_cubic_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_catenary>(csound, "hc_funicular_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_diocles_curve>(csound, "hc_cissoid_curve", "i", "i", csnd::thread::i); // alias for diocles
    csnd::plugin<cs_gaussian>(csound, "hc_gaussian_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_typed>(csound, "hc_transeg_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_mouth>(csound, "hc_kiss_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_bicorn>(csound, "hc_cocked_hat_curve", "i", "i", csnd::thread::i);

    // Core
    csnd::plugin<cs_segment>(csound, "hc_segment", "i", "iii", csnd::thread::i);
    csnd::plugin<cs_gen>(csound, "hc_hypercurve", "i", "iiim", csnd::thread::i);
    csnd::plugin<cs_gen>(csound, "hc_gen", "i", "iiim", csnd::thread::i);

    //deprecated
    csnd::plugin<cs_normalize_curve>(csound, "hc_normalize_y", "", "iii", csnd::thread::i);
}
