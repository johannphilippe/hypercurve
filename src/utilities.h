/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include<cmath>
#include<iostream>
#include<vector>
#include<memory>
#include<cstring>
#include<functional>
#include"fpng/src/fpng.h"

namespace hypercurve {

inline int round(double x)
{
    return (x - std::floor(x)) ? 1 : 0;
}

inline double pos(double x)
{
    return (x > 0) ? x : 0;
}

inline double frac(double x1, double x2)
{
    return double(x1) / double(x2);
}

inline double limit(double min, double max, double v)
{
    if(v > max) return max;
    if(v < min) return min;
    return v;
}

template<typename T>
std::shared_ptr<T> share(T t)
{
    return std::make_shared<T>(t);
}

struct point
{
    point(double x_, double y_)
        : x(x_)
        , y(y_)
    {}

    point() {}

    template<typename T>
    double distance_to(T &p)
    {
        const double a2 = (p.x - x) * (p.x - x);
        const double b2 = (p.y - y) * (p.y - y);
        const double dist = sqrt(a2 + b2);
        return dist;
    }

    void print()
    {
        printf("point \n\tx = %f \n\ty = %f\n", x, y);
    }

    double x,y;
};
using control_point = point;

struct curve_point : public point
{
    curve_point(float x_, float y_, float curve_ = 0.0f) :
        point(x_, y_), curve(curve_)
    {}

    curve_point(point &p, float curve_ = 0.0f) :
        point(p), curve(curve_)
    {}

   float curve;
};

// x must be between 0 and 1, 0 being the y1 x position, 1 the y2 x position
inline double linear_interpolation(double y1, double y2, double x)
{
    return y1 + (x * (y2 - y1)) ;
}

// x must be between 0 and 1, 0 being the y1 x position, 1 the y2 x position
// Not tested
inline double cubic_interpolation(double y1, double y2, double x)
{
    return y1 + (std::pow(x,3) * (y2 - y1));
}

inline double cubic_root(double x) {return std::pow(x, 1.0 / 3.0);}
inline double squared(double x) {return x*x;}
inline double cubed(double x) {return x*x*x;}

// Linear or exp interpolation, based on Csound GEN16 algorithm
template<typename T>
inline T log_exp_point(T beg, T ending, int dur, int idx, double typ)
{
    if(typ == 0) {
        return beg + (idx * (( ending - beg) / dur));
    }else {
        double type = limit(-10, 10, typ);
        return beg + (ending - beg) * (1 - std::exp(idx * type / (dur - 1))) / (1 - std::exp(type));
    }
}

// returns x between 0 and 1 (1 being x == x2, 0 being x == x1)
// make sure x1 <= x <= x2
inline double relative_position(double x1, double x2, double x)
{
    if(!(x >= x1 && x <= x2) || !(x1 < x2))
        throw(std::runtime_error("Make sure x1 <= x <= x2 and x1 < x2"));
    const double factor = 1.0 / (x2 - x1);
    return (x * factor) - (x1 * factor);
}

std::vector<double> linspace(size_t size)
{
    std::vector<double> v(size);
    for(size_t i = 0; i < size; i++)
        v[i] = (double(i) / double(size));
    return v;
}

inline double pi_angle(double angle)
{
    return angle * M_PI / 180.0;
}

// Pi relative angle
//https://www.desmos.com/calculator/aqmmojldge?lang=fr
inline double rotate(double x, double y, double angle, std::function<double(double)> f)
{
    return (f(x * std::cos(angle) - y * std::sin(angle) ) - x * std::sin(angle) ) / std::cos(angle);
}

// Windows
inline double hanning(int index, int length) {
    return  0.5 * (1 - cos(2 * M_PI * index / (length - 1 )));
}

inline double hamming(int index, int length) {
    return 0.54 - (0.46 * cos(2 * M_PI * index / (length - 1)));
}

inline double blackman(int index, int length) {
    return 0.42 - (0.5 * cos(2 * M_PI * index / (length - 1))) + (0.08 * cos(4 * M_PI * index / (length - 1)));
}

// To allow external allocators
template<typename T>
struct memory_vector
{
   struct iterator : public std::iterator<std::forward_iterator_tag, T>
   {
      iterator() : ptr(nullptr) {}
      iterator(T *ptr_)
         : ptr(ptr_)
      {}
      iterator(const iterator&o)
         : ptr(o.ptr)
      {}
      ~iterator() {
      }
      iterator& operator=(const iterator&o) {ptr = o.ptr; return *this;}
      iterator& operator++() {++ptr; return *this;} //prefix increment
      iterator& operator++(int) {return ptr++;}
      T * operator->() {return ptr;}
      T& operator*() const {return *ptr;}
      bool operator==(const iterator& i) {return ptr == i.ptr;}
      bool operator!=(const iterator& i) {return ptr != i.ptr;}

      T *ptr;
   };

   memory_vector() {}
   memory_vector(const memory_vector<T>& other)
   {
       _size = other._size;
       _data = new T[_size];
       std::copy(other._data, other._data + other._size, _data );
   }
   memory_vector(size_t size_)
   {
      _size = size_;
      _data = new T[_size];
   }
   memory_vector(T *ptr, size_t size_)
   {
      _size = size_;
      _data = ptr;
   }

   // Copy std::vector
   memory_vector(std::vector<T> &v)
   {
      _size = v.size();
      _data = new T[_size];
      std::copy(v.begin(), v.end(), _data);
   }

   // Move std::vector - allow allocation with :
   // memory_vector<type>({arg1, arg2, arg3...});
   memory_vector(std::vector<T> &&v)
   {
      _size = v.size();
      _data = new T[_size];
      std::move(v.begin(), v.end(), _data);
   }

   // operators

   memory_vector<T> operator=(const memory_vector<T> &other)
   {
       return memory_vector<T>(other);
   }

   void init(T *ptr, size_t size_)
   {
       _size = size_;
       _data = ptr;
   }

    ~memory_vector()
   {
       delete[] _data;
   }

   void resize(size_t size_)
   {
      if(_size > 0)
         delete[] _data;
      _size=  size_;
      _data = new T[_size];
   }

   T& operator[](size_t index)
   {
      return _data[index];
   }

   size_t size() {return _size;}

   iterator begin() {return iterator(_data);}
   iterator end() {return iterator(_data + _size);}
   T* data() {return _data;}

   T *_data;
   size_t _size;
};


// PNG utils

//using color = std::array<uint8_t, 4>;
struct color : public std::array<uint8_t, 4>
{
};

constexpr const color white{{255, 255, 255, 255}};
constexpr const color black{{0, 0, 0, 255}};
constexpr const color red{{192, 9,9, 255}};
constexpr const color purple{{106, 9, 192, 255}};


struct png
{
    png(size_t width_ = 2048, size_t height_ = 1024, color background_ = black, color foreground_ = purple)
        : width(width_)
        , height(height_)
        , background(background_)
        , foreground(foreground_)
        , data(width * height, background)
    {}

    void set(size_t x, size_t y, color c)
    {
        data[ width * (height - y - 1) + x] = c;
    }

    void set_curve_point(double x, double y)
    {
        size_t ix = x * width;
        size_t iy = y * height;
        set(ix, iy, foreground);
    }

    void fill_curve_point(double x, double y, bool waveform = false)
    {
        size_t ix = x * width;
        size_t iy = y * height; // y position of point

        size_t half = height / 2;
        if(!waveform)
        {
            for(size_t i = 0; i < iy; ++i)
            {
                set(ix, i, foreground);
            }
        } else
        {
            for(size_t i = std::min(iy, half); i < std::max(iy, half); ++i)
            {
                set(ix, i, foreground);
            }
        }
    }

    void draw_grid(size_t xdiv = 10, size_t ydiv = 10, color c = white)
    {
        for(size_t i = 1; i < xdiv; ++i)
        {
            size_t yb = 0;
            size_t ye = height;
            size_t x = i * frac(width, xdiv);
            for(size_t y = yb; y < ye; ++y)
            {
                set(x, y, c);
            }
        }
        for(size_t i = 1; i < ydiv; ++i)
        {
            size_t xb = 0;
            size_t xe = width;
            size_t y = i * frac(height, ydiv);
            for(size_t x = xb; x < xe; ++x)
            {
                set(x, y, c);
            }
        }

    }

    void draw_curve(double *samples, size_t size, bool fill = true, bool waveform = false )
    {
        for(size_t i = 0; i < size; ++i)
        {
            double x = frac(i,size);
            double y = (waveform) ? (samples[i] + 1.0) / 2.0 : samples[i];
            if(fill)
                fill_curve_point(x, y, waveform);
            else
                set_curve_point(x, y);
        }

    }

    void write_as_png(std::string path)
    {
        fpng::fpng_encode_image_to_file(path.c_str(), (void *) data.data(), width, height, 4);
    }

protected:
    size_t width, height;
    color background, foreground;
    std::vector<color> data;
};


}

#endif // UTILITIES_H
