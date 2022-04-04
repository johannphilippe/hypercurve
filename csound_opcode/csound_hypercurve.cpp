#include"../src/hypercurve.h"
//#include"../src/curve_lib.h"
#include<plugin.h>
#include<OpcodeBase.hpp>
#include<unordered_map>

using namespace hypercurve;
/////////////////////////////////////////////////////////////
//
// Csound Hypercurve Plugin
//
// * Thoughts
// If we implement the full thing as opcodes, for each curve it
// would require an opcode that fits the particular curve.
// On the other side a JSON implementation could be good to
// ease the process, but wouldn't allow to write it in Csound
// (So no live coding of curves)
//
/////////////////////////////////////////////////////////////

void *mtx;


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

protected:
    size_t _definition;
    double y_start;
    memory_vector<double> samples;
    memory_vector<double>::iterator it;
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
      std::cout << "curve index : " << index << std::endl;
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

    int index;

    std::shared_ptr<control_point> cp1, cp2;

    std::shared_ptr<catmull_rom_spline_curve> _curve;
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
      std::cout << "segment index : " << index << std::endl;
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
      mtx = csound->get_csound()->Create_Mutex(0);
      csound::LockGuard g(csound->get_csound(), mtx);
      _mem.allocate(csound, (int)inargs[0]);
      _initialize(inargs[0], inargs[1], _mem.data());
      process_init();
      for(size_t i = 2; i < in_count(); ++i)
      {
          if(segment_map.find(int(inargs[i])) == segment_map.end())
              throw(std::runtime_error("Segment not found in data"));

          y_start = (double) ((i > 2) ?  segment_map[inargs[i - 1] ]->y_destination : inargs[1]);
          segment_map[inargs[i]]->set_y_start(y_start);
          process_one(segment_map[inargs[i]]);
      }

      index = curve_map.map(this);

      // Weirdly, it works below
      AsciiPlotter p("helloplot" , 80, 15);
      p.addPlot(std::vector<double>(linspace(inargs[0])), std::vector<double>(this->samples.begin(), this->samples.end()), "helloplot", '*');
      p.show();

      outargs[0] = index;
      return OK;
  }

  int deinit()
  {
      //curve_map.unmap(index);
      return OK;
  }

  int index;
  double y_start;
  csnd::AuxMem<double> _mem;
};

// Do the same as below with trigger (trig_run_hypercurve) with frequency and trigger parameter.
struct cs_run_curve : public csnd::Plugin<1, 2>
{
  int init()
  {
      index = inargs[0];
      if(curve_map.find(index) == curve_map.end())
          throw(std::runtime_error("This curve does not exist"));
      return OK;
  }

  int kperf()
  {
      phasor = inargs[1];
      i_phasor = floor(phasor * curve_map[index]->get_definition());
      outargs[0] = limit(-1, 1, curve_map[index]->get_sample_at(i_phasor));
      return OK;
  }

  int aperf()
  {
      for(size_t i = 0; i < ksmps(); ++i)
      {
          phasor = inargs.vector_data<MYFLT>(0)[i];
          i_phasor = floor(phasor * curve_map[index]->get_definition());
          outargs(0)[i] = limit(-1, 1, curve_map[index]->get_sample_at(i_phasor));

      }

      return OK;
  }

private:
  double phasor;
  size_t i_phasor;
  int index;
};


#include <modload.h>

void csnd::on_load(Csound *csound) {
    std::cout << "loading csound hypercurve" << std::endl;
    // Helpers
    csnd::plugin<cs_control_point>(csound, "control_point", "i", "ii", csnd::thread::i);

    // Curve types
    csnd::plugin<cs_diocles_curve>(csound, "diocles_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_cubic_curve>(csound, "cubic_curve", "i", "", csnd::thread::i);
    csnd::plugin<cs_quad_bezier>(csound, "quadratic_bezier_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_cub_bezier>(csound, "cubic_bezier_curve", "i", "ii", csnd::thread::i);
    csnd::plugin<cs_catmull_rom>(csound, "catmull_rom_curve", "i" , "ii", csnd::thread::i);

    // Core
    csnd::plugin<cs_segment>(csound, "segment", "i", "iii", csnd::thread::i);
    csnd::plugin<cs_curve>(csound, "hypercurve", "i", "iim", csnd::thread::i);

    // Run
    csnd::plugin<cs_run_curve>(csound, "run_hypercurve", "k", "ik", csnd::thread::ik);
    csnd::plugin<cs_run_curve>(csound, "run_hypercurve", "a", "ia", csnd::thread::ia);
}
