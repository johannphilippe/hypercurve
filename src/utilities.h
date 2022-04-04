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

namespace hypercurve {

double pos(double x)
{
    return (x > 0) ? x : 0;
}

double frac(double x1, double x2)
{
    return double(x1) / double(x2);
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

// returns x between 0 and 1 (1 being x == x2, 0 being x == x1)
// make sure x1 <= x <= x2
inline double relative_position(double x1, double x2, double x)
{
    if(!(x >= x1 && x <= x2) || !(x1 < x2))
        throw(std::runtime_error("Make sure x1 <= x <= x2 and x1 < x2"));
    const double factor = 1.0 / (x2 - x1);
    return (x * factor) - (x1 * factor);
}

inline double limit(double min, double max, double v)
{
    if(v > max) return max;
    if(v < min) return min;
    return v;
}

std::vector<double> linspace(size_t size)
{
    std::vector<double> v(size);
    for(size_t i = 0; i < size; i++)
        v[i] = (double(i) / double(size));
    return v;
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
      for(size_t i = 0; i < v.size(); ++i)
          _data[i] = v[i];
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
}

#endif // UTILITIES_H
