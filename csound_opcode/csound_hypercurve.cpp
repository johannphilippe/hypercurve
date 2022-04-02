#include"../src/hypercurve.h"
#include<plugin.h>
#include<unordered_map>
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


using namespace hypercurve;

template<typename T>
struct increment_map
{
    increment_map()
    {
    }
    int map(T to_map)
    {
        _map[_index] = to_map;
        ++_index;
        return _index - 1;
    }

    T& operator[](int idx)
    {
        return _map[idx];
    }

    int _index = 0;
    T _map[100];
};

static increment_map< std::shared_ptr<curve_base> > curve_base_map;
static increment_map< segment > segment_map;
static increment_map< curve > curve_map;

struct cs_diocles_curve : public csnd::Plugin<1, 1>
{
  int init()
  {
      _curve.allocate(csound, 1);
      std::cout << "init " << std::endl;
      _curve[0] = diocles_curve(inargs[0]);
      std::cout << "diocles ok  " << std::endl;
      int index = curve_base_map.map(share(_curve[0]));
      printf("index mapped \n");
      outargs[0] = index;
      printf("retur ok \n");
      return OK;
  }

private:
  csnd::AuxMem<diocles_curve> _curve;
};

struct cs_segment : public csnd::Plugin<3, 1>
{
  int init()
  {
      index = segment_map.map(segment(inargs[0], inargs[1], curve_base_map[inargs[2]]));
      outargs[0] = index;
      return OK;
  }

private:
  int index;
};

struct cs_curve : public csnd::Plugin<3, 1>
{
  int init()
  {
      size_t len = inargs.vector_data<MYFLT>(2).len();
      _segments.allocate(csound, len);
      for(size_t i = 0; i < len; ++i)
          _segments[i] = segment_map[i];
      index = curve_map.map(curve(inargs[0], inargs[1], _segments.data(), len));
      outargs[0] = index;
      return OK;
  }
  int index;
  csnd::AuxMem<segment> _segments;
};

struct cs_run_curve : public csnd::Plugin<2, 1>
{
  int init()
  {
      index = inargs[0];
      return OK;
  };

  int kperf()
  {
      double phasor = inargs[1];
      size_t i_phasor = floor(phasor * curve_map[index].get_definition());
      outargs[0] = curve_map[index].get_sample_at(i_phasor);

      return OK;
  };
private:
  int index;
};


struct simple_thing : csnd::Plugin<1,0>
{

    int init()
    {
        outargs[0]= (rand() % 126) / double(126.0);
        return OK;
    }
};

#include <modload.h>

void csnd::on_load(Csound *csound) {
    std::cout << "loading csound hypercurve "<< std::endl;
    csnd::plugin<cs_diocles_curve>(csound, "diocles_curve", "i", "i", csnd::thread::i);
    csnd::plugin<cs_segment>(csound, "segment", "i", "iii", csnd::thread::i);
    csnd::plugin<cs_curve>(csound, "hypercurve", "i", "iii", csnd::thread::i);
    csnd::plugin<cs_run_curve>(csound, "run_hypercurve", "k", "ik", csnd::thread::ik);
}
