/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef HYPERCURVE_FAUST_H
#define HYPERCURVE_FAUST_H

/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef HYPERCURVE_H
#define HYPERCURVE_H

/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef CORE_H
#define CORE_H

#include<cmath>
#include<vector>
#include<memory>
/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef CURVE_LIB_H
#define CURVE_LIB_H

/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include<cmath>
#include<iostream>
#include<vector>
#include<unordered_map>
#include<memory>
#include<cstring>
#include<functional>
// fpng.h - unlicense (see end of fpng.cpp)

#include <stdlib.h>
#include <stdint.h>
#include <vector>

#ifndef FPNG_TRAIN_HUFFMAN_TABLES
	// Set to 1 when using the -t (training) option in fpng_test to generate new opaque/alpha Huffman tables for the single pass encoder.
	#define FPNG_TRAIN_HUFFMAN_TABLES (0)
#endif

namespace fpng
{
	// ---- Library initialization - call once to identify if the processor supports SSE.
	// Otherwise you'll only get scalar fallbacks.
	void fpng_init();

	// ---- Useful Utilities

	// Returns true if the CPU supports SSE 4.1, and SSE support wasn't disabled by setting FPNG_NO_SSE=1.
	// fpng_init() must have been called first, or it'll assert and return false.
	bool fpng_cpu_supports_sse41();

	// Fast CRC-32 SSE4.1+pclmul or a scalar fallback (slice by 4)
	const uint32_t FPNG_CRC32_INIT = 0;
	uint32_t fpng_crc32(const void* pData, size_t size, uint32_t prev_crc32 = FPNG_CRC32_INIT);

	// Fast Adler32 SSE4.1 Adler-32 with a scalar fallback.
	const uint32_t FPNG_ADLER32_INIT = 1;
	uint32_t fpng_adler32(const void* pData, size_t size, uint32_t adler = FPNG_ADLER32_INIT);

	// ---- Compression
	enum
	{
		// Enables computing custom Huffman tables for each file, instead of using the custom global tables. 
		// Results in roughly 6% smaller files on average, but compression is around 40% slower.
		FPNG_ENCODE_SLOWER = 1, 
		
		// Only use raw Deflate blocks (no compression at all). Intended for testing.
		FPNG_FORCE_UNCOMPRESSED = 2,
	};

	// Fast PNG encoding. The resulting file can be decoded either using a standard PNG decoder or by the fpng_decode_memory() function below.
	// pImage: pointer to RGB or RGBA image pixels, R first in memory, B/A last.
	// w/h - image dimensions. Image's row pitch in bytes must is w*num_chans.
	// num_chans must be 3 or 4. 
	bool fpng_encode_image_to_memory(const void* pImage, uint32_t w, uint32_t h, uint32_t num_chans, std::vector<uint8_t>& out_buf, uint32_t flags = 0);

#ifndef FPNG_NO_STDIO
	// Fast PNG encoding to the specified file.
	bool fpng_encode_image_to_file(const char* pFilename, const void* pImage, uint32_t w, uint32_t h, uint32_t num_chans, uint32_t flags = 0);
#endif

	// ---- Decompression
		
	enum
	{
		FPNG_DECODE_SUCCESS = 0,				// file is a valid PNG file and written by FPNG and the decode succeeded
		
		FPNG_DECODE_NOT_FPNG,					// file is a valid PNG file, but it wasn't written by FPNG so you should try decoding it with a general purpose PNG decoder

		FPNG_DECODE_INVALID_ARG,				// invalid function parameter

		FPNG_DECODE_FAILED_NOT_PNG,				// file cannot be a PNG file
		FPNG_DECODE_FAILED_HEADER_CRC32,		// a chunk CRC32 check failed, file is likely corrupted or not PNG
		FPNG_DECODE_FAILED_INVALID_DIMENSIONS,  // invalid image dimensions in IHDR chunk (0 or too large)
		FPNG_DECODE_FAILED_DIMENSIONS_TOO_LARGE, // decoding the file fully into memory would likely require too much memory (only on 32bpp builds)
		FPNG_DECODE_FAILED_CHUNK_PARSING,		// failed while parsing the chunk headers, or file is corrupted
		FPNG_DECODE_FAILED_INVALID_IDAT,		// IDAT data length is too small and cannot be valid, file is either corrupted or it's a bug

		// fpng_decode_file() specific errors
		FPNG_DECODE_FILE_OPEN_FAILED,
		FPNG_DECODE_FILE_TOO_LARGE,
		FPNG_DECODE_FILE_READ_FAILED,
		FPNG_DECODE_FILE_SEEK_FAILED
	};

	// Fast PNG decoding of files ONLY created by fpng_encode_image_to_memory() or fpng_encode_image_to_file().
	// If fpng_get_info() or fpng_decode_memory() returns FPNG_DECODE_NOT_FPNG, you should decode the PNG by falling back to a general purpose decoder.
	//
	// fpng_get_info() parses the PNG header and iterates through all chunks to determine if it's a file written by FPNG, but does not decompress the actual image data so it's relatively fast.
	// 
	// pImage, image_size: Pointer to PNG image data and its size
	// width, height: output image's dimensions
	// channels_in_file: will be 3 or 4
	// 
	// Returns FPNG_DECODE_SUCCESS on success, otherwise one of the failure codes above.
	// If FPNG_DECODE_NOT_FPNG is returned, you must decompress the file with a general purpose PNG decoder.
	// If another error occurs, the file is likely corrupted or invalid, but you can still try to decompress the file with another decoder (which will likely fail).
	int fpng_get_info(const void* pImage, uint32_t image_size, uint32_t& width, uint32_t& height, uint32_t& channels_in_file);

	// fpng_decode_memory() decompresses 24/32bpp PNG files ONLY encoded by this module.
	// If the image was written by FPNG, it will decompress the image data, otherwise it will return FPNG_DECODE_NOT_FPNG in which case you should fall back to a general purpose PNG decoder (lodepng, stb_image, libpng, etc.)
	//
	// pImage, image_size: Pointer to PNG image data and its size
	// out: Output 24/32bpp image buffer
	// width, height: output image's dimensions
	// channels_in_file: will be 3 or 4
	// desired_channels: must be 3 or 4 
	// 
	// If the image is 24bpp and 32bpp is requested, the alpha values will be set to 0xFF. 
	// If the image is 32bpp and 24bpp is requested, the alpha values will be discarded.
	// 
	// Returns FPNG_DECODE_SUCCESS on success, otherwise one of the failure codes above.
	// If FPNG_DECODE_NOT_FPNG is returned, you must decompress the file with a general purpose PNG decoder.
	// If another error occurs, the file is likely corrupted or invalid, but you can still try to decompress the file with another decoder (which will likely fail).
	int fpng_decode_memory(const void* pImage, uint32_t image_size, std::vector<uint8_t>& out, uint32_t& width, uint32_t& height, uint32_t& channels_in_file, uint32_t desired_channels);

#ifndef FPNG_NO_STDIO
	int fpng_decode_file(const char* pFilename, std::vector<uint8_t>& out, uint32_t& width, uint32_t& height, uint32_t& channels_in_file, uint32_t desired_channels);
#endif

	// ---- Internal API used for Huffman table training purposes

#if FPNG_TRAIN_HUFFMAN_TABLES
	const uint32_t HUFF_COUNTS_SIZE = 288;
	extern uint64_t g_huff_counts[HUFF_COUNTS_SIZE];
	bool create_dynamic_block_prefix(uint64_t* pFreq, uint32_t num_chans, std::vector<uint8_t>& prefix, uint64_t& bit_buf, int& bit_buf_size, uint32_t *pCodes, uint8_t *pCodesizes);
#endif

} // namespace fpng

#include <random>
#include<string>
#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif

namespace hypercurve {

template<typename Numeric, typename Generator = std::mt19937>
Numeric random(Numeric from, Numeric to)
{
    thread_local static Generator gen(std::random_device{}());

    using dist_type = typename std::conditional
    <
        std::is_integral<Numeric>::value
        , std::uniform_int_distribution<Numeric>
        , std::uniform_real_distribution<Numeric>
    >::type;

    thread_local static dist_type dist;

    return dist(gen, typename dist_type::param_type{from, to});
}

inline int round(double x)
{
    return (x - std::floor(x)) ? 1 : 0;
}

inline double pos(double x)
{
    return (x > 0) ? x : 0;
}

template<typename T, typename TT>
inline double fraction(T x1, TT x2)
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
inline std::shared_ptr<T> share(T t)
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

    point operator/(point &other) {return point{this->x/other.x, this->y/other.y};}
    point operator*(point &other) {return point{this->x*other.x, this->y*other.y};}
    point operator+(point &other) {return point{this->x+other.x, this->y+other.y};}
    point operator-(point &other) {return point{this->x-other.x, this->y-other.y};}
    point &operator--() {--x; --y; return *this;}
    point &operator--(int) {--x; --y; return *this;}
    point &operator++() {++x; ++y; return *this;}
    point &operator++(int) {++x; ++y; return *this;}

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
        throw(std::runtime_error( std::string("Make sure x1 <= x <= x2 and x1 < x2 \nx1 = ")
                                  + std::to_string(x1) + " & x = "
                                  + std::to_string(x) + " & x2 = "
                                  + std::to_string(x2)));
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

constexpr static const double hamming_scaling_constant = 0.08;
constexpr static const double hamming_scaling_factor = 1.0 / (1.0 - hamming_scaling_constant);
inline double hamming(int index, int length) {
    return 0.54 - (0.46 * cos(2 * M_PI * double(index) / ( double(length) - 1.)));
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
      iterator operator+(int add) {return ptr+add;}
      iterator operator-(int add) {return ptr-add;}
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

    ~memory_vector()
   {
       delete[] _data;
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

   T max() {
       T mx = 0;
       for(size_t i = 0; i < size(); ++i)
           if(_data[i] > mx) mx = _data[i];
       return mx;
   }
   T min() {
       T mx = 0;
       for(size_t i = 0; i < size(); ++i)
           if(_data[i] < mx) mx = _data[i];
       return mx;
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
struct color //: public std::array<uint8_t, 4>
{
    uint8_t r,g,b,a;
};

constexpr const color white{255, 255, 255, 255};
constexpr const color black{0, 0, 0, 255};
constexpr const color red{192, 9,9, 255};
constexpr const color purple{106, 9, 192, 255};

struct png
{
    png(size_t width_ = 2048, size_t height_ = 1024,
        color background_ = black, color foreground_ = purple)
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
            size_t x = i * hypercurve::fraction(width, xdiv);
            for(size_t y = yb; y < ye; ++y)
            {
                set(x, y, c);
            }
        }
        for(size_t i = 1; i < ydiv; ++i)
        {
            size_t xb = 0;
            size_t xe = width;
            size_t y = i * hypercurve::fraction(height, ydiv);
            for(size_t x = xb; x < xe; ++x)
            {
                set(x, y, c);
            }
        }

    }

    void draw_curve(double *samples, size_t size,
                    bool fill = true, bool waveform = false )
    {
        for(size_t i = 0; i < size; ++i)
        {
            double x = hypercurve::fraction(i,size);
            double y = (waveform) ? (samples[i] + 1.0) / 2.0 : samples[i];
            if(fill)
                fill_curve_point(x, y, waveform);
            else
                set_curve_point(x, y);
        }

    }

    void write_as_png(std::string path)
    {
        fpng::fpng_encode_image_to_file(path.c_str(),
                                        (void *) data.data(), width, height, 4);
    }

protected:
    size_t width, height;
    color background, foreground;
    std::vector<color> data;
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

    bool has(int index)
    {
        return this->find(index) != this->end();
    }

    void dump()
    {
        std::cout << "index now is : " << _index << std::endl;
        for(size_t i = 1; i < _index; ++i)
        {
            std::cout << "incr map : " << i << " ptr : " << &this->at(i) << std::endl;
        }
    }

    size_t _index = 1;
};

}

#endif // UTILITIES_H

#include<cmath>
#include<iostream>
#include<functional>
#include<vector>
/*=============================================================================
   Copyright (c) 2021 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef SPLINE_H
#define SPLINE_H

#include<iostream>
#include<vector>
#include<unordered_map>
#include<cmath>
#include<functional>
#include<memory.h>

using namespace std;
namespace hypercurve {

template<typename T>
class cubic_spline
{
public:

    cubic_spline()
    {
    }

    ~cubic_spline()
    {
    }

    void gauss_elimination_ls(int m, int n)
    {
        int i, j, k;
        for( i = 0; i < m - 1; i++)
        {
            for(k = i + 1; k < m; k++)
            {
                double term = tridiagonal[k][i] / tridiagonal[i][i];
                for(j = 0; j < n; j++)
                {
                    tridiagonal[k][j] = tridiagonal[k][j] -
                            term * tridiagonal[i][j];
                }
            }
        }

        // Back substitution
        for(i = m - 1; i >= 0; i--)
        {
            sig_temp[i] = tridiagonal[i][n - 1];
            for(j = i + 1; j < n - 1; j++)
            {
                sig_temp[i] = sig_temp[i] - tridiagonal[i][j] * sig_temp[j];
            }
            sig_temp[i] = sig_temp[i] / tridiagonal[i][i];
        }
    }

    void cs_coeff_calculation(int n)
    {
        for(int i = 0; i < n; i++)
        {
            d[i] = y[i];
            b[i] = sig[i] / 2.0;
            a[i] = (sig[i + 1] - sig[i]) / (h[i] * 6.0);
            c[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * sig[i] + sig[i + 1]) / 6.0;

        }
    }

    void tridiagonal_cubic_splin_gen(int n)
    {
        int i;
        for(i = 0; i < n - 1; i++)
        {
            tridiagonal[i][i] = 2 * (h[i]+ h[i + 1]);
        }
        for(i = 0; i < n - 2; i++)
        {
            tridiagonal[i][i + 1] = h[i + 1];
            tridiagonal[i + 1][i] = h[i + 1];
        }
        for(i = 1; i < n; i++)
        {
            tridiagonal[i - 1][n - 1] = (y[i + 1] - y[i]) * 6.0 / h[i] -
                    (y[i] - y[i - 1]) * 6.0 / h[i - 1];
        }
    }

    T equation(T xval, T x_i, int i)
    {
        return (a[i] * pow((xval - x_i), 3)) + (b[i] * pow((xval - x_i), 2)) +
                (c[i] * (xval - x_i)) + d[i];
    }

    void calculate_interpolation(int)
    {
        int key_index = 0;
        T key = x[key_index];
        for (int i = 0; i < interp.size(); i++)
        {
            T xval = T(i) / T(n_precision);
            if((xval < x[0]) || (xval > x.back()) ) continue;
            if(xval > x[key_index + 1]) {
                key_index++;
                key = x[key_index];
            }
            interp[i] = equation(xval, key, key_index);  //equation_map[key](xval);
        }
       // equation_map.clear();
    }

    void resize(int n)
    {
        x.resize(n + 1);
        y.resize(n + 1);

        h.resize(n);

        a.resize(n);
        b.resize(n);
        c.resize(n);
        d.resize(n);

        sig.resize(n + 1);
        sig_temp.resize(n - 1);

        tridiagonal.resize(n - 1);
        for(int i = 0; i < tridiagonal.size(); i++)
        {
            tridiagonal[i].resize(n);
        }
    }

    void reset(int n, vector<point> p)
    {
        for(int i = 0; i < p.size(); i++)
        {
            x[i] = p[i].x;
            y[i] = p[i].y;
        }

        ::memset(h.data(), 0, n * sizeof(T));
        ::memset(a.data(), 0, n * sizeof(T));
        ::memset(b.data(), 0, n * sizeof(T));
        ::memset(c.data(), 0, n * sizeof(T));
        ::memset(d.data(), 0, n * sizeof(T));
        ::memset(sig.data(), 0, (n + 1) * sizeof(T)  );
        ::memset(sig_temp.data(), 0, (n - 1) * sizeof(T) );

        for(int i = 0; i < tridiagonal.size(); i++)
        {
            ::memset(tridiagonal[i].data(), 0, n * sizeof(T) );
        }

        ::memset(interp.data(), 0, sizeof(T) * n_precision);
    }

    void interpolate(int n)
    {

        int i;
        for(i = 0; i < n; i++)
            h[i] = x[i + 1] - x[i];

        sig[0] = 0;
        sig[n] = 0;

        tridiagonal_cubic_splin_gen(n);
        gauss_elimination_ls(n -1, n);

        for(i = 1; i < n; i++)
            sig[i] = sig_temp[i - 1];

        cs_coeff_calculation(n);
        calculate_interpolation(n);
    }

    vector<T>& interpolate_from_points(vector<point> &p,
                                            int precision, point)
    {
        int n = p.size() - 1;
        if(n_points != p.size())
        {
            this->resize(n);
            n_points = p.size();
        }
        if(precision != n_precision)
        {
            interp.resize(precision);
            n_precision = precision;
        }

        reset(n, p);

        this->interpolate(n);

        return interp;
    }

protected:
    int n_points = 0, n_precision = 0;
    vector<T> x, y, interp, a, b, c, d, sig, sig_temp, h;
    vector< vector<T> > tridiagonal;
};

}

#endif // SPLINE_H

namespace  hypercurve {

//////////////////////////////////////////////////
// Base Curve
//////////////////////////////////////////////////
class curve_base
{
public:
    curve_base() {}
    virtual ~curve_base() {}

    // This method is the one to implement/override for simple curves (who only need x and constants to calculate y)
    inline virtual double process(double x) {return x;};
    // When you want to retrieve a single sample from a curve, it is recommanded to use this one. It is not necessary to override it for simple curves,
    // but may be necessary for complex (bezier, spline, catmullrom ...)
    inline virtual double process(size_t i, size_t size) {return process(hypercurve::fraction(i, size));}
    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it)
    {
        memory_vector<double>::iterator begin_ptr = it;
        double max = 0.0;
        for(size_t i = 0; i < size; ++i)
        {
            const double x = hypercurve::fraction(i, size);
            double res = scale(x);
            if( std::abs(res) > max) max = std::abs(res);
            //if(inverted) res = process_invert(x, res);
            *it = res;
            ++it;
        }

        if(inverted)
        {
            memory_vector<double> tmp(definition);
            for(size_t i = 0; i < definition; ++i)
            {
                // Trouver la valeur en X correspondant à un y linéaire donné

                point diff(1, y_destination - y_start);
                double x_itp = fraction(i, definition);
                point perp_origin(0 + x_itp, diff.y);
                point perp_dest(1 + x_itp, -diff.y);
                bool found = false;
                size_t start_x = i;
                double y = linear_interpolation(perp_dest.y, perp_origin.y, fraction(i, definition));
                while(!found)
                {
                    start_x --;
                    double y_itp = linear_interpolation(perp_dest.y, perp_origin.y, fraction(start_x, definition) );

                    // Faux
                    double curve_y = *(begin_ptr + start_x);

                    std::cout << y << " " << y_itp << " " << curve_y << std::endl;

                    if(( y_itp > curve_y && y < curve_y) || (y_itp < curve_y && y > curve_y))
                    {
                        found = true;
                        tmp[i] = curve_y;
                        break;

                    }
                    y = y_itp;
                }

            }

            for(size_t i = 0; i < definition; ++i)
                *(begin_ptr+i) = tmp[i];
        }

        return max;
    }

    // Do not override this (or make sure to implement members definition)
    inline virtual void init(double y_start_, double y_dest_, size_t definition_) {
        y_start = y_start_;
        y_destination = y_dest_;
        definition = definition_;
        abs_diff = std::abs(y_start - y_destination);
        offset = std::min(y_start, y_destination);
        on_init();
    };

    // Override this one insted of init to avoid y_start and y_destination
    // affectation repetition
    inline virtual void on_init() {}

    // Allows inversion of curve (y symetry on a linear x_start/x_end axis).
    bool inverted = false;
protected:

    inline virtual double scale(double x)
    {
        if(y_start > y_destination)
            return process(1.0 - x) * abs_diff + offset;
        return process(x) * abs_diff + offset;
    }

    inline double process_invert(double x, double y)
    {

        const double lin = linear_interpolation(y_start, y_destination, x);
        return lin + (lin-y); //return lin + (lin - y) ;
    }

    size_t definition;
    double y_start, y_destination, abs_diff, offset;
};

using linear_curve = curve_base;

// Allows you to make a symetry on a x_start/y_start x_end/y_end linear axis
inline std::shared_ptr<curve_base> invert(std::shared_ptr<curve_base> cb)
{
    cb->inverted = true;
    return cb;
}

//////////////////////////////////////////////////
// Curve Library
//////////////////////////////////////////////////

class diocles_curve : public curve_base
{
public:
    diocles_curve(double a_)
        : a(a_)
        , compensation(1. / process_diocles(1.0) )
    {}
    inline double process(double x) override
    {
        return process_diocles(x) * (compensation);
    }

private:

    inline double process_diocles(double x)
    {
         return std::sqrt( std::pow(x, 3) / (2 * a - x) );
    }

    double a, compensation;
};

using cissoid_curve = diocles_curve;

class cubic_curve: public curve_base
{
public:
    inline double process(double x) override
    {
        return std::pow(x, 3);
    }
};

// Choose the exponent of X
class power_curve  : public curve_base
{
public:
    power_curve(double exponent_)
        : exponent(exponent_)
    {}
    inline double process(double x) override
    {
        return std::pow(x, exponent) ;
    }

private:
    double exponent;
};

class hanning_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        return  hanning(x * definition, definition * 2);
    }
};

class hamming_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        return  (hamming(x * definition, definition * 2) - hamming_scaling_constant) * hamming_scaling_factor;
    }
};

class blackman_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        return  blackman(x * definition, definition * 2);
    }
};

// See https://mathcurve.com/courbes2d/bouche/bouche.shtml
class mouse_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        const double _x = x -1;
        return std::sqrt( std::pow( (a*a) - (_x*_x), 3 ) / std::pow(a, 4) );
    }

    constexpr static const double a = 1.0;
};
using kiss_curve = mouse_curve;

// See https://mathcurve.com/courbes2d/bicorne/bicorne.shtml
class bicorn_curve : public curve_base
{
public:
    // Sign true = positive, false = negative
    bicorn_curve(bool sign_)
        : sign(sign_ == true ? 1 : -1)
    {}

    inline double process(double x) override
    {
        const double _x = x-1;
        if(sign == 1)
            return process_bicorn(_x) / process_bicorn(0);
        return process_bicorn(_x);
    }

    constexpr static const double a = 1.0;
    int sign;

protected:
    inline double process_bicorn(double x)
    {
        return  ((a*a) - (x*x)) / (2*a + (std::sqrt(a*a - x*x) * sign));
    }
};
using cocked_hat_curve = bicorn_curve;

//////////////////////////////////////////////////
// Typed curve - inspired by Csound GEN 16
//////////////////////////////////////////////////
class typed_curve : public curve_base
{
public:
    typed_curve(double type_)
        : type(-type_)
    {}

    inline double process(double x) override
    {
        return log_exp_point<double>(0, 1, definition, x * definition, type);
    }

    double type;
};

//////////////////////////////////////////////////
// Gauss Curve
//////////////////////////////////////////////////
class gauss_curve :  public curve_base
{
public:
    gauss_curve(const double A_, const double c_)
        : A(A_)
        , c(c_)
    {
        process_width();
    }

    inline double process(double x) override
    {
        const double _x = (x * half_width) - half_width;
        double gauss = A * std::exp(-std::pow(_x, 2) / (2 * (c*c)));
        gauss -= y_offset;
        gauss /= (A - y_offset); //  scale between offset and 1
        return gauss;
    }

protected:

    void process_width()
    {
        half_width = std::sqrt(2 * std::log(10) * c);
        y_offset = A * std::exp(-std::pow(half_width, 2) / (2 * (c*c)));
    }

    double A, c;
    double half_width, y_offset;

};
using gaussian_curve = gauss_curve;

//////////////////////////////////////////////////
// Catenary (funicular) curve
// In french "courbe de chainette"
// https://mathcurve.com/courbes2d.gb/chainette/chainette.shtml
// The more 'a' is important, the more this curve is straight
//////////////////////////////////////////////////
class catenary_curve : public curve_base
{
public:
    catenary_curve(double a_)
        : a(a_)
        , height(catenary_process(1))
    {
        if(!(a>0))
            throw (std::runtime_error("Catenary curve : a must be > 0"));
    }
    double process(double x) override
    {
        return catenary_process(x) / height;
    }
protected:
    double catenary_process(double x)
    {
        return a * (std::cosh(x/a)) - a;
    }

    double a, height;

    constexpr static const double mult = 12;
};
using funicular_curve = catenary_curve;

//////////////////////////////////////////////////
// Conchal curve
// https://mathcurve.com/courbes2d.gb/conchale/conchale.shtml
// TODO invert x y axis
// x must be scaled from -a to max_x
// Must be scaled and axis inverted
// Still to do
// Angle b = -4.713
//////////////////////////////////////////////////
class conchal_curve : public curve_base
{
public:
    conchal_curve(double a_, double c_)
        : a(a_)
        , c(c_)
        , w(std::sqrt((a*a) + (c*c)))
        , height(a + w)
    {
        if( (a > c) || (a <= 0) || (c <= 0) || (c <  (a * 1.2)))
            throw(std::runtime_error("Conchal curve : a must be < c, c must be >= a * 1.2, a and c must be > 0"));
        max_x = a + (std::sqrt((a*a) + (c*c)));
    }

    double process(double x) override
    {
        const double step = height / double(this->definition);
        const double mx = x * height - a + step;
        const double my = conchal_process(mx) / conchal_process(-a);
        std::cout << "mx : " << mx << " my : " << my << std::endl;
        const double res = (rotate(mx, my, pi_angle(-90), [&](double x_) {
           return this->conchal_process(x_);
        }) + a ) / (height);

        //std::cout << "res = " << res << std::endl;
        if(std::isnan(res))
            std::cout <<"NAN ALERT at " << x  << "  mx : " << mx << std::endl;
        return res;
    }

protected:
    inline double conchal_process(double x)
    {
        std::cout <<"conchal : " << ((c*c) + (a*a) - (x*x) )
                    * ((c*c) - (a*a) + (x*x)) << "   " << (x+a) << std::endl;
        return std::sqrt(
                    ((c*c) + (a*a) - (x*x))
                    * ((c*c) - (a*a) + (x*x))
                    ) / (x+a);
    }

    inline double scaled_conchal(double x)
    {
        double _x = (x * (max_x + a)) - a;
        double _conchal_res = conchal_process(_x);
        return x + (_conchal_res + a);
    }

    double a, c, max_x, w, height;
};

//////////////////////////////////////////////////
// Tightrope walker cure
// https://mathcurve.com/courbes2d.gb/danseur/danseur.shtml
// The more a is sperior to b, the more is is "straight"
// abs(b) is the curve x max point
// Also a must not be too close to b (undefined behavior) : a = 1.01 b = 1 is ok. a = 1.0001 is not.
//////////////////////////////////////////////////

class tightrope_walker_curve : public curve_base
{
public:
    tightrope_walker_curve(double a_, double b_)
        : a(a_)
        , b( std::abs(b_) )
    {
        if( (a < 0) || (a <= b) )
            throw(std::runtime_error("a must be > abs(b), a must be > 0"));
    }

    void on_init() override
    {
        max_x = b - b / double(definition);
        max_y = abs(process_tightrope_walker(max_x));
    }

    double process(double x) override
    {
        if( ( x * b) >= max_x) return process_tightrope_walker(max_x) / max_y;
        return process_tightrope_walker(x * b) / max_y;
    }

protected:
    inline double process_tightrope_walker(double x)
    {
        return (x * (a-x)) / std::sqrt((b*b)-(x*x));
    }

    double a,b, max_x, max_y;
};

//////////////////////////////////////////////////
// Toxoid curve
// AKA duplicatrix cubic
// https://mathcurve.com/courbes2d/cubicduplicatrice/cubicduplicatrice.shtml
// a needs to be <= 0 (see constructor)
//////////////////////////////////////////////////

class toxoid_curve : public curve_base
{

public:
    toxoid_curve(double a_)
        : a(-std::abs(a_))
        , height(toxoid_process( 1))
    {
    }

    double process(double x) override
    {
        return toxoid_process(x) / height;
    }

protected:

    inline double toxoid_process(double x)
    {
        return std::sqrt(x * std::pow(x - a/2, 2));
        //return std::sqrt( ( (x*x) * (x - a/2) ) / a  );
    }

    double a, height;
};
using duplicatrix_cubic = toxoid_curve;

//////////////////////////////////////////////////
// user defined curves
//
// callback should return an y value between 0 and 1
// for each x between 0 and 1
//////////////////////////////////////////////////
class user_defined_curve : public curve_base
{
public:
    user_defined_curve() {}
    user_defined_curve(std::function<double(double)> f)
        : callback(f)
    {}

    inline double process(double x) override
    {
        return callback(x);
    }
protected:
    std::function<double(double)> callback;
};

//////////////////////////////////////////////////
// vararg polynomial
//
// If you give three parameters a, b, c
// it will return ax^3 + bx^2 + cx
//////////////////////////////////////////////////
class polynomial_curve : public curve_base
{
public:
    polynomial_curve(memory_vector<double> args)
        : constants(args)
    {}

    polynomial_curve(std::vector<double> args)
        : constants(args)
    {}

    polynomial_curve(double *args, size_t size)
        : constants(args, size)
    {}

    double process(double x) override
    {
        double res = 0;
        for(size_t i = 0; i < constants.size(); ++i)
        {
            res += std::pow(x, constants.size() - i) * constants[i];
        }
        // then scale.
        return res;
    }

    memory_vector<double> constants;
};

//////////////////////////////////////////////////
// Bezier Curves
//////////////////////////////////////////////////

class bezier_curve_base : public curve_base
{
public:
    inline double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        double max = 0.0;
        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = process_bezier(double(cnt) / double(size));
        r2 = process_bezier(double(cnt + 1) / double(size));

        for(size_t i = 0; i < size ; i++)
        {
            const double x = double(i) / double(size);
            while(cnt < size)
            {
                if( (r1.first <= x ) && (r2.first >= x))
                    break;
                r1 = process_bezier( double(cnt) / double(size) );
                r2 = process_bezier( double(cnt + 1) / double(size) );
                ++cnt;
            }

            double relative_x = relative_position(r1.first, r2.first, x);
            double linear_interp = linear_interpolation(r1.second, r2.second, relative_x);
            if(std::abs(linear_interp) > max) max = std::abs(linear_interp);
            if(inverted) linear_interp = process_invert(x, linear_interp);
            *it = linear_interp;
            ++it;
        }

        return max;
    }

    inline virtual double process(double) override
    {
        throw(std::runtime_error("Unimplemented for Bezier curve"));
    }

    // This one should be implemented instead of the above one (if wanted to process single bezier point)
    inline virtual double process(size_t i, size_t size) override
    {
        return process(hypercurve::fraction(i, size));
    }
protected:
    inline virtual std::pair<double, double> process_bezier(double x) {return {x,0};}
};

//https://math.stackexchange.com/questions/3280087/can-the-t-of-a-quadratic-bezier-curve-be-found-given-a-x-or-y-using
class quadratic_bezier_curve : public bezier_curve_base
{
public:
    quadratic_bezier_curve(control_point cp)
        : _control_point(cp)
    {}

    void on_init() override
    {
        c_x = (3 * (_control_point.x));
        b_x = (-c_x);
        a_x = (1 - c_x - b_x);
        c_y = ( 3 * (_control_point.y - y_start));
        b_y = (-c_y);
        a_y = (y_destination - y_start - c_y - b_y);
    }

private:
    inline std::pair<double, double> process_bezier(double x) override
    {
        return {
            (a_x * std::pow(x, 3)) + (b_x * std::pow(x, 2)) + (c_x * x),
            (a_y * std::pow(x,3)) + (b_y * std::pow(x,2)) + (c_y * x) + y_start
        };
    }

    control_point _control_point;
    double c_x, b_x, a_x, c_y, b_y, a_y;
};

// To get direct interpolation for Bezier
//https://stackoverflow.com/questions/15505392/implementing-ease-in-update-loop/15506642#15506642
// t is :
// 0.25 = 0.6*t(1-t)^2 + 0.5*t^2(1-t) + t^3
// 0.25 = 0.6 * t *(1 - t) ^ 2 + 0.5 * t^2 * (1 - t) + t ^ 3

// Or
// https://stackoverflow-com.translate.goog/questions/5883264/interpolating-values-between-interval-interpolation-as-per-bezier-curve?_x_tr_sl=en&_x_tr_tl=fr&_x_tr_hl=fr&_x_tr_pto=sc

class cubic_bezier_curve : public bezier_curve_base
{
public:
    cubic_bezier_curve(control_point _cp1, control_point _cp2)
        : _control_point1(_cp1)
        , _control_point2(_cp2)
    {}

    void on_init() override
    {
        c_x = ((-1) * 0 + 3 * _control_point1.x - 3 * _control_point2.x + 1);
        b_x = (3 * 0 - 6 * _control_point1.x + 3 * _control_point2.x);
        a_x = ((-3) * 0 + 3 * _control_point1.x);
        c_y = ((-1) * y_start + 3 * _control_point1.y - 3 * _control_point2.y + y_destination);
        b_y = (3 * y_start - 6 * _control_point1.y + 3 * _control_point2.y);
        a_y = ((-3) * y_start + 3 * _control_point1.y);
    }

private:
    std::pair<double, double> process_bezier(double x) override
    {
        return {
            (std::pow(x, 3) * c_x)
                    + (std::pow(x, 2) * b_x)
                    + (x * a_x),
            (std::pow(x, 3) * c_y)
                    + (std::pow(x, 2) * b_y)
                    + (x * a_y)
                    + y_start
        };
    }
    control_point _control_point1, _control_point2;
    double c_x, b_x, a_x, c_y, b_y, a_y;
};

//////////////////////////////////////////////////
// Spline Curves
//////////////////////////////////////////////////

class cubic_spline_curve : public virtual curve_base
{
public:
    cubic_spline_curve(std::vector<point> cp)
        : _control_points( std::move(cp) )
    {
        if(_control_points.size() < 3)
            throw(std::runtime_error("Control point list size must be >= 3"));
    }

    void on_init() override
    {
       _control_points.insert(_control_points.begin(), curve_point(0, y_start));
        // For each point, determine absolute position from relative position
        for(size_t i = 1; i < _control_points.size(); i++)
        {
            _control_points[i].x += _control_points[i-1].x;
        }
    }

    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        std::vector<double>& res = spl.interpolate_from_points(_control_points, size, point{1.0, 1.0});
        double max = 0.0;
        for(size_t i = 0; i < size; i++)
        {
            if( std::abs(res[i]) > max ) max = std::abs(res[i]);
            if(inverted) res[i] = process_invert(hypercurve::fraction(i, size), res[i]);
            *it = res[i];
            ++it;
        }
        return max;
    }

private:
    cubic_spline<double> spl;
    std::vector<point> _control_points;
};

// User passes control points P0 and P3 assuming that P1(0,y) and P2(1,y). Calculation will be relative, and will be  rescaled after.
// Based on https://www.desmos.com/calculator/552cpvzfxw?lang=fr

// To get a real centripetal or chordal, we should be based on https://www.desmos.com/calculator/9kazaxavsf?lang=fr
// alpha =  Parametric constant: 0.5 for the centripetal spline, 0.0 for the uniform spline, 1.0 for the chordal spline.
class catmull_rom_spline_curve : public curve_base
{
public:
    catmull_rom_spline_curve(point p0, point p3)
        : _cp0(p0)
        , _cp3(p3)
    {}

    void on_init() override
    {
        _cp1 = point(0, y_start);
        _cp2 = point(1, y_destination);
    }

    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = process_catmul_rom(0);
        r2 = process_catmul_rom(1.0 / double(size));
        double max = 0.0;
        for(size_t i =0; i < size; i++)
        {
            const double x = double(i) / double(size);
            while(cnt < size)
            {
                if((r1.first <= x) && (r2.first >= x))
                    break;
                r1 = process_catmul_rom( double(cnt) / double(size) );
                r2 = process_catmul_rom( double(cnt + 1) / double(size) );
                ++cnt;
            }

            double relative_x = relative_position(r1.first, r2.first, x);
            double linear_interp = linear_interpolation(r1.second, r2.second, relative_x);
            if(std::abs(linear_interp) > max) max = std::abs(linear_interp);
            if(inverted) linear_interp = process_invert(x, linear_interp);
            *it = linear_interp;
            ++it;
        }
        return  max;
    }
private:

    std::pair<double, double> process_catmul_rom(double x)
    {
        const double rx = 0.5 * ((_cp1.x * 2.0) + (-_cp0.x + _cp2.x) * x
                                 + ((_cp0.x * 2.0) - (_cp1.x * 5.0)
                                    + (_cp2.x * 4.0) - _cp3.x ) * (x*x)
                                 + (-_cp0.x  + (_cp1.x * 3.0) - (_cp2.x * 3.0) +_cp3.x) * (x*x*x));
        const double ry = 0.5 * ((_cp1.y * 2.0) + (-_cp0.y + _cp2.y) * x
                                 + ((_cp0.y * 2.0) - (_cp1.y * 5.0)
                                    + (_cp2.y * 4.0) - _cp3.y ) * (x*x)
                                 + (-_cp0.y  + (_cp1.y * 3.0) - (_cp2.y * 3.0) +_cp3.y) * (x*x*x));
        return {rx, ry};
    }

    control_point _cp0, _cp3, _cp1, _cp2;
};

class lagrange_polynomial_curve : public curve_base
{

public:
    lagrange_polynomial_curve(std::vector<control_point> pts)
        : c_pts(pts.size() + 2)
    {
        std::copy(pts.begin(), pts.end(), c_pts.data() + 1);
    }
    lagrange_polynomial_curve(memory_vector<control_point> pts)
        : c_pts(pts.size() + 2)
    {
        std::copy(pts.begin(), pts.end(), c_pts.data() + 1);
    }

    // Used for preallocated situations (e.g. Csound)
    lagrange_polynomial_curve(control_point *pts, size_t size)
        : c_pts(pts, size)
    {}

    void on_init() override
    {
        c_pts[0] = control_point(0, y_start);
        c_pts[c_pts.size() - 1] = control_point(1, y_destination);
    }

    inline double process(double x) override
    {
        return process_lagrange(x);
    }

protected:
    double process_lagrange(double x)
    {
        yp = 0;
        for(size_t i = 0; i < c_pts.size(); ++i)
        {
            p = 1;
            for(size_t j = 0; j < c_pts.size(); ++j)
            {
                if(j != i)
                    p = p * (x - c_pts[j].x)/( c_pts[i].x  - c_pts[j].x);
            }
            yp = yp + p * c_pts[i].y;
        }
        return yp;
    }

    double p, yp;
    memory_vector<control_point> c_pts;
};

}
#endif // CURVE_LIB_H

/*
MIT License

Copyright (c) 2017  Joe Hood

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include<vector>
#include<string>

#define MAX_CURVES 10

using namespace std;

class AsciiPlotter
{
private:
	static const char _default_markers[10];
	int _width;
	int _height;
	int _curves;
	string _title;
	string _xlabel;
	string _ylabel;
	bool _legend;
	vector<double> _ydata[MAX_CURVES];
	string _labels[MAX_CURVES];
	char _markers[MAX_CURVES];
public:
	AsciiPlotter();
	AsciiPlotter(string title);
	AsciiPlotter(string title, int width, int height);
	~AsciiPlotter();
	void plot(const char *plotfile, const char *datafile);
    void addPlot( vector<double> ydata, string label, char marker);
    void show();
	void xlabel(string label);
	void ylabel(string label);
	void legend();
};

#include<iostream>
#include<utility>
namespace hypercurve {
///////////////////////////////////////////////////
// The segment class
////////////////////////////////////////////////////

class segment
{
public:
    segment(double frac, double y_dest, std::shared_ptr<curve_base> c)
        : fractional_size(frac)
        , y_destination(y_dest)
        , _curve( std::move(c) )
    {}

    segment() {}

    virtual double process(memory_vector<double>::iterator &it, size_t size)
    {
        return _curve->process_all(size, it);
    }

    virtual void init(double y_start_, size_t definition)
    {
        y_start = y_start_;
        _curve->init(y_start, y_destination, definition);
    }
    void rescale_x(double factor)
    {
        fractional_size *= factor;
    }

    double fractional_size;
    double y_destination;
    double y_start = 0;
    std::shared_ptr<curve_base> _curve;
protected:
};

/////////////////////////////////////////////////////
// Curve class
////////////////////////////////////////////////////

class curve
{
public:
    curve(size_t definition_, double y_start_, std::vector< segment > segs_)
        : definition(definition_)
        , y_start(y_start_)
        , segs( std::move(segs_) )
        , samples(definition)
    {
        check_total_size();
        init();
        process();
    }

    curve() {}
    virtual  ~curve() {}

    virtual void init()
    {
        segs[0].init(y_start, std::round(segs[0].fractional_size * definition) );
        for(size_t i = 1; i < segs.size(); i++)
        {
            segs[i].init(segs[i - 1].y_destination, std::round(segs[i].fractional_size * definition) );
        }
    }

    void process()
    {
       memory_vector<double>::iterator it = samples.begin();
        for(size_t i = 0; i < segs.size(); i++)
        {
            // For each segment, we must give it a real size (size_t), and an iterator position
            size_t seg_size = std::round(segs[i].fractional_size * definition);
            double seg_max = segs[i].process(it, seg_size);
            if(std::abs(seg_max)  > max) max = std::abs(seg_max);
        }
        find_extremeness();
    }

    void normalize_y(double target_min, double target_max)
    {
        find_extremeness();

        for(size_t i = 0; i < samples.size(); i++)
        {
            samples[i] = ((samples[i] - min ) / ambitus )  * std::abs(target_max - target_min) + target_min;
        }
    }

    std::pair<double, double> find_extremeness()
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
        return std::make_pair(min, max);
    }

    void ascii_display(std::string name, std::string label, char marker)
    {
        AsciiPlotter plot(name, 80, 15);
        plot.addPlot(std::vector<double>(samples.data(), samples.data() + definition), label, marker);
        plot.legend();
        plot.show();
    }

    double *get_samples() {return samples.data();}
    double get_sample_at(size_t i) {return samples[i];}
    size_t get_definition() {return definition;}

    // Operators

    double& operator[](size_t index)  { return samples[index];}

    curve& operator *=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] *= other.get_sample_at(i);
        return *this;
    }
    curve &operator *=(double k)
    {
        for(auto & it : samples)
            it *= k;
        return *this;
    }

    curve& operator /=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] /= other.get_sample_at(i);
        return *this;
    }
    curve &operator /=(double k)
    {
        for(auto & it : samples)
            it /= k;
        return *this;
    }

    curve operator *(curve &c2)
    {
        curve c(*this);
        c *= c2;
        return c;
    }

    curve operator *(double k)
    {
        curve c(*this);
        c *= k;
        return c;
    }

    curve operator /(curve &c2)
    {
        curve c(*this);
        c /= c2;
        return c;
    }

    curve operator / (double k)
    {
        curve c(*this);
        c /= k;
        return c;
    }

    curve& operator+=(double k)
    {
        for(auto & it : samples)
            it += k;
        return *this;
    }

    curve& operator +=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] += other.get_sample_at(i);
        return *this;
    }

    curve operator +(double k)
    {
        curve c(*this);
        c += k;
        return c;
    }
    curve operator +(curve &other)
    {
        curve c(*this);
        c += other;
        return c;
    }

    curve& operator -=(double k)
    {
        for(auto & it : samples)
            it -= k;
        return *this;
    }
    curve& operator -=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] -= other.get_sample_at(i);
        return *this;
    }

    curve operator -(double k)
    {
        curve c(*this);
        c -= k;
        return c;
    }

    curve operator -(curve &other)
    {
        curve c(*this);
        c -= other;
        return c;
    }

protected:
    void check_total_size()
    {
        double x = 0;
        for(auto & it : segs)
            x += it.fractional_size;
        if( x != 1.0 )
        {
            this->rescale(x);
        }
    }

    void rescale(double x)
    {
        double factor = (1. / x);
        for(auto & it : segs) {
            it.rescale_x(factor);
        }
    }

    double max = 0.0, min = 0.0, ambitus = 0.0;
    size_t definition;
    double y_start;
    memory_vector< segment > segs;
    memory_vector<double> samples;
};

}
#endif // CORE_H

/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef MODULATOR_LIB_H
#define MODULATOR_LIB_H

#include<cmath>
#include<variant>
#include<type_traits>

namespace hypercurve {
////////////////////////////////////////////////////
// Interpolator is a purely algorithmic curve
// it is used to scale modulators
////////////////////////////////////////////////////
class interpolator
{
public:
    interpolator(double y_start, std::vector<point> pts)
        : itps(std::move(pts))
    {
        // add point 0
        itps.insert(itps.begin(), point(0, y_start));
        for(size_t i = 1; i < itps.size(); ++i)
            itps[i].x += itps[i-1].x;
    }

    virtual double interpolate(double x)
    {
        int seg = which_segment(x);
        double relative_x = relative_position(itps[seg].x, itps[seg + 1].x, x);
        return process(itps[seg].y, itps[seg + 1].y, relative_x);
    }

    virtual double process(double y1, double y2, double x)
    {
        return linear_interpolation(y1, y2, x);
    }

private:
    int which_segment(double x)
    {
        for(size_t i = 0; i < itps.size() - 1; ++i)
        {
            if(x >= itps[i].x && x <  itps[i + 1].x)
                return i;
        }
        return -1;
    }
    std::vector<point> itps;
};

using linear_interpolator = interpolator;

class cubic_interpolator : public interpolator
{
public:
    cubic_interpolator(double y_start ,std::vector<point> pts)
        : interpolator( y_start, std::move(pts) )
    {}
    double process(double y1, double y2, double x) override
    {
        return cubic_interpolation(y1, y2, x);
    }
};

////////////////////////////////////////////////////
// Amplitude can be static (double between 0 and 1)
// or dynamic (interpolator)
////////////////////////////////////////////////////

class amplitude {
public:
    virtual double get_amplitude(double x) {return x;}
};

class amplitude_fixed : amplitude
{
public:
    amplitude_fixed(double d)
        : _amplitude(d)
    {}
    double get_amplitude(double) override {return _amplitude;}
    double _amplitude;
};

class amplitude_interpolated : amplitude
{
public:
    amplitude_interpolated(std::shared_ptr<interpolator> itp)
        : _amplitude(std::move(itp))
    {}
    double get_amplitude(double x) override {return _amplitude->interpolate(x);}
    std::shared_ptr<interpolator> _amplitude;
};

////////////////////////////////////////////
// Modulators are like curve_base
// They are special kind of "curve" that take
// an amplitude parameter (e.g. noise, oscillator)
///////////////////////////////////////////

template<typename Amp>
class modulator_base : public curve_base, public Amp
{
public:
    modulator_base(double amp)
        : Amp(amp)
    {}

    modulator_base(shared_ptr<interpolator> itp)
        : Amp(itp)
    {}

    double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        double max = 0.0;
        for(size_t i = 0; i < size; ++i)
        {
            const double x = hypercurve::fraction(i, size);
            double res = Amp::get_amplitude(x) * process(x);
            if( std::abs(res) > max) max = std::abs(res);
            *it = res;
            ++it;
        }
        return max;
    }

protected:
};

template<typename Amp = amplitude_fixed>
class noise_modulator : public modulator_base<Amp>
{
public:
    noise_modulator(double amp, size_t precision = 256)
        : modulator_base<Amp>(amp)
        , _precision(precision)
    {}

    noise_modulator(shared_ptr<interpolator> itp, size_t precision = 256)
        : modulator_base<Amp>(itp)
        , _precision(precision)
    {}

    inline double process(double) override
    {
        return ( (rand() % (_precision*2) ) - _precision) / double(_precision);
    }

protected:
    size_t _precision;
};

template<typename Amp>
class sine_modulator : public modulator_base<Amp>
{
public:
    sine_modulator(double amp, double freq)
        : modulator_base<Amp>(amp)
        , _freq(freq)
    {}
    sine_modulator(std::shared_ptr<interpolator> itp, double freq)
        : modulator_base<Amp>(itp)
        , _freq(freq)
    {}

    inline double process(double x) override
    {
        return std::sin(x * M_PI * 2 * _freq);
    }
private:
    double _freq;

};

//////////////////////////////////////////////////
// Chebyshev
// T=1 is stable (basically scaled between -1 and 1)
// T=2 is not, and shouldn't be used if you don't
// know what you're doing
//////////////////////////////////////////////////

template<typename Amp, int T = 1>
class chebyshev_modulator : public modulator_base<Amp>
{
public:
    chebyshev_modulator(double amp, int n_)
        : modulator_base<Amp> (amp)
        , n(n_)
    {}
    chebyshev_modulator(shared_ptr<interpolator> tip, int n_)
        : modulator_base<Amp> (tip)
        , n(n_)
    {}

    inline double process(double x) override
    {

        const double t = std::acos(( x * 2.0) - 1.0 );
        if constexpr(T == 1)
        {
            return std::cos(n * t);
        } else // T == 2
        {
            return std::sin( (n + 1) * t) / std::sin(t);
        }
    }
private:

    const double n;
};
}

#endif // MODULATOR_LIB_H

/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef EXTRA_H
#define EXTRA_H

#include<iostream>
#include<unordered_map>

#include<memory>

namespace hypercurve {

//////////////////////////////////////////////////////////////////////
// Random and algorithmic curve composition
//////////////////////////////////////////////////////////////////////

// A map of curve types
enum curve_base_index {
    linear_i = 0,
    diocles_i = 1,
    cubic_i = 2,
    power_i = 3,
    hanning_i = 4,
    hamming_i = 5,
    blackman_i = 6,
    gaussian_i = 7,
    toxoid_i = 8,
    catenary_i = 9,
    tightrope_walker_i = 10,
    cubic_bezier_i = 11,
    quadratic_bezier_i = 12,
    cubic_spline_i = 13,
    catmull_rom_spline_i = 14,
    polynomial_i = 15,
    user_defined_i = 16,
    typed_i = 17,
    mouse_i = 18,
    bicorn_i = 19,
    lagrange_polynomial_i = 20,

    // keep that last to get size
    size_i
};

const char *get_curve_base_name (const curve_base_index b)
{
    switch(b) {
    case linear_i: return "linear";
    case diocles_i: return "diocles";
    case cubic_i: return "cubic";
    case power_i: return "power";
    case hanning_i: return "hanning";
    case hamming_i: return "hamming";
    case blackman_i: return "blackman";
    case gaussian_i: return "gaussian";
    case toxoid_i: return "toxoid";
    case catenary_i: return "catenary";
    case tightrope_walker_i: return "tightrope_walker";
    case cubic_bezier_i: return "cubic_bezier";
    case quadratic_bezier_i: return "quadratic_bezier";
    case cubic_spline_i: return "cubic_spline";
    case catmull_rom_spline_i: return "catmull_rom_spline";
    case polynomial_i: return "polynomial";
    case user_defined_i: return "user_defined";
    case typed_i: return "typed";
    case mouse_i: return "mouse";
    case bicorn_i: return "bicorn";
    case lagrange_polynomial_i: return "lagrange_polynomial";
    default: return "";
    }
    return "";
}

std::shared_ptr<curve_base> get_curve_from_index(curve_base_index n,
                                                 std::vector<double> args,
                                                 std::vector<control_point> cps)
{
    switch (n) {
    case linear_i: return share(linear_curve());
    case diocles_i: return share(diocles_curve(args[0]));
    case cubic_i: return share(cubic_curve());
    case power_i: return share(power_curve(args[0]));
    case hanning_i: return share(hanning_curve());
    case hamming_i: return share(hamming_curve());
    case blackman_i: return share(blackman_curve());
    case gaussian_i: return share(gaussian_curve(args[0], args[1]));
    case toxoid_i: return share(toxoid_curve(args[0]));
    case catenary_i: return share(catenary_curve(args[0]));
    case tightrope_walker_i: return share(tightrope_walker_curve(args[0], args[1]));
    case cubic_bezier_i: return share(cubic_bezier_curve(cps[0], cps[1]));
    case quadratic_bezier_i: return share(quadratic_bezier_curve(cps[0]));
    case catmull_rom_spline_i: return share(catmull_rom_spline_curve(cps[0], cps[1]));
    case polynomial_i: return share(polynomial_curve(args));
    case cubic_spline_i: return share(cubic_spline_curve(cps));
    case typed_i: return share(typed_curve(args[0]));
    case mouse_i: return share(mouse_curve());
    case bicorn_i: return share(bicorn_curve(args[0] > 0 ));
    case lagrange_polynomial_i: return share(lagrange_polynomial_curve(cps));
    default: return share(linear_curve());
    }
}

std::pair<std::vector<double>, std::vector<control_point>> random_args_generator(curve_base_index n)
{
    std::vector<double> args;
    std::vector<control_point> cps;

    switch (n) {
    case linear_i: return {{},{}};
    case diocles_i: return {{random<double>(0.50001, 10)}, {}};
    case cubic_i: return {{},{}};
    case power_i: return {{(double)random<int>(1, 32)}, {}};
    case hanning_i: return {{},{}};
    case hamming_i: return {{},{}};
    case blackman_i: return {{},{}};
    case gaussian_i: return {{random<double>(0.1, 4), random<double>(0.00001, 4)}, {}};
    case toxoid_i: return {{random<double>(0.0001, 10)}, {}};
    case catenary_i: return {{random<double>(0.0001, 1.9999)},{}};
    case tightrope_walker_i: {
        double a = random<double>(1, 3);
        double b = a - random<double>(0.01, 0.99);
        return {{a, b}, {}};

    };
    case cubic_bezier_i: {
        std::vector<control_point> cps{
          control_point(random<double>(0, 1), random<double>(0, 1)),
          control_point(random<double>(0, 1), random<double>(0, 1))
        };
        std::sort(cps.begin(),cps.end(), [&](control_point f1,control_point f2) {
            return f1.x < f2.x;
        });

        return {{}, cps};
    };
    case quadratic_bezier_i: return {{}, {control_point(random<double>(0, 1),
                                                        random<double>(0, 1))}};
    case catmull_rom_spline_i: {
        std::vector<control_point> cps {
          control_point(random<double>(0.01, 3) * -1.0, random<double>(0.01, 34) * -1),
          control_point(random<double>(1.01, 3) , random<double>(1.01, 3) ),
        };
        return {{}, cps};
    };
    case polynomial_i: {
        size_t nargs = (rand() % 10) + 1;
        std::vector<double> args(nargs);
        for(size_t i = 0; i < nargs; ++i)
            args[i] = random<double>(0, 10);
        return {args, {}};
        };
    case lagrange_polynomial_i:
    case cubic_spline_i: {
        size_t nargs = (rand() % 4) + 3;
        std::vector<control_point> cps(nargs);
        for(size_t i = 0; i < nargs; ++i)
            cps[i] = control_point(random<double>(0, 1), random<double>(0, 1));
        std::sort(cps.begin(), cps.end(), [&](control_point cp1, control_point cp2) {
           return cp1.x < cp2.x;
        });
        return {{}, cps};
    };
    case typed_i: return {{random<double>(-10, 10)}, {}};
    case mouse_i: return {{},{}};
    case bicorn_i: return {{(double)random<int>(-1, 1)}, {}};

    default: return {{},{}};
    }
}

std::pair<curve, std::string> random_curve_composer( size_t max_segs = 16, int min = 0, int max = 1,
                             size_t definition = 4096, bool envelop = false, bool waveform = false,
                             bool force_curve_type = false, curve_base_index forced = linear_i)
{
    auto gen_curve = [&](){
        return force_curve_type ? forced : static_cast<curve_base_index>(rand() % static_cast<int>(size_i));
    };

    std::string cnames = "";
    size_t nsegs = 1 + (rand() % max_segs);
    if(nsegs < 2) nsegs = 2;

    std::vector<segment> segs(nsegs);
    for(size_t i = 0; i < nsegs; ++i) {
        double frac_size = random<double>(0.1, 1);

        double dest = waveform ? random<double>(-1, 1) : random<double>(0, 1);
        curve_base_index index =  gen_curve();

        while( index == user_defined_i || (envelop && ((index == polynomial_i) || (index == lagrange_polynomial_i) || (index == cubic_spline_i))))
            index = gen_curve();

        cnames += std::string(get_curve_base_name(index)) + "_";
        std::pair<std::vector<double>, std::vector<control_point>> args = random_args_generator(index);
        if(i == (nsegs - 1) && (envelop | waveform) )
            dest = 0;
        segs[i] = segment(frac_size, dest, get_curve_from_index(index, args.first, args.second) );
    }

    double y_start = (envelop | waveform) ? 0 : random<double>(0, 1);
    curve c(definition, y_start, segs);
    c.normalize_y(min, max);
    return {c, cnames};
}

}

#endif // EXTRA_H

#endif // HYPERCURVE_H

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
static int hc_normalize_y(int crv, double min, double max)
{
    faust_curve_map[crv]->normalize_y(min, max);
    return crv;
}

// Just inverts a curve base on its linear axis (symetry on linear axis)
static int hc_invert_curve_base(int cb)
{
    invert(faust_curve_base_map[cb]);
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

// Curve bases functions

static const int max_segments = 64;

static int hc_diocles_curve(double a) {return faust_curve_base_map.map(share(diocles_curve(a)));}
static int hc_linear_curve(int) {return faust_curve_base_map.map(share(linear_curve()));}
static int hc_cubic_curve(int) {return faust_curve_base_map.map(share(cubic_curve()));}
static int hc_power_curve(double exponent) {return faust_curve_base_map.map(share(power_curve(exponent)));}
static int hc_hanning_curve(int) {return faust_curve_base_map.map(share(hanning_curve()));}
static int hc_hamming_curve(int) {return faust_curve_base_map.map(share(hamming_curve()));}
static int hc_blackman_curve(int) {return faust_curve_base_map.map(share(blackman_curve()));}
static int hc_gaussian_curve(double A, double c) {return faust_curve_base_map.map(share(gaussian_curve(A,c)));}
static int hc_toxoid_curve(double a) {return faust_curve_base_map.map(share(toxoid_curve(a)));}
static int hc_catenary_curve(double a) {return faust_curve_base_map.map(share(catenary_curve(a)));}
static int hc_tightrope_walker_curve(double a, double b) {return faust_curve_base_map.map(share(tightrope_walker_curve(a, b)));}

static int hc_typed_curve(double t) {return faust_curve_base_map.map(share(typed_curve(t)));}
static int hc_mouse_curve(int) {return faust_curve_base_map.map(share(mouse_curve()));}
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
// fpng.cpp 1.0.6 - Fast 24/32bpp .PNG image writer/reader. See unlicense at the end of this file.
// PNG's generated by this code have been tested to load successfully with stb_image.h, lodepng.cpp, wuffs, libpng, and pngcheck.
//
// Uses code from the simple PNG writer function by Alex Evans, 2011. Released into the public domain: https://gist.github.com/908299
// Some low-level Deflate/Huffman functions derived from the original 2011 Google Code version of miniz (public domain by R. Geldreich, Jr.): https://code.google.com/archive/p/miniz/
// Low-level Huffman code size function: public domain, originally written by: Alistair Moffat, alistair@cs.mu.oz.au, Jyrki Katajainen, jyrki@diku.dk, November 1996.
//
// Optional config macros:
// FPNG_NO_SSE - Set to 1 to completely disable SSE usage, even on x86/x64. By default, on x86/x64 it's enabled.
// FPNG_DISABLE_DECODE_CRC32_CHECKS - Set to 1 to disable PNG chunk CRC-32 tests, for improved fuzzing. Defaults to 0.
// FPNG_USE_UNALIGNED_LOADS - Set to 1 to indicate it's OK to read/write unaligned 32-bit/64-bit values. Defaults to 0, unless x86/x64.
//
// With gcc/clang on x86, compile with -msse4.1 -mpclmul -fno-strict-aliasing
// Only tested with -fno-strict-aliasing (which the Linux kernel uses, and MSVC's default).
//

#include <assert.h>
#include <string.h>

#ifdef _MSC_VER
	#pragma warning (disable:4127) // conditional expression is constant
#endif

// Set FPNG_NO_SSE to 1 to completely disable SSE usage.
#ifndef FPNG_NO_SSE
	#define FPNG_NO_SSE (0)
#endif

// Detect if we're compiling on x86/x64
#if defined(_M_IX86) || defined(_M_X64) || defined(__i386__) || defined(__i386) || defined(__i486__) || defined(__i486) || defined(i386) || defined(__ia64__) || defined(__x86_64__)
	#define FPNG_X86_OR_X64_CPU (1)
#else
	#define FPNG_X86_OR_X64_CPU (0)
#endif

#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE
	#ifdef _MSC_VER
		#include <intrin.h>
	#endif
	#include <xmmintrin.h>		// SSE
	#include <emmintrin.h>		// SSE2
	#include <smmintrin.h>		// SSE4.1
	#include <wmmintrin.h>		// pclmul
#endif

#ifndef FPNG_NO_STDIO
	#include <stdio.h>
#endif

// Allow the disabling of the chunk data CRC32 checks, for fuzz testing of the decoder
#ifndef FPNG_DISABLE_DECODE_CRC32_CHECKS
	#define FPNG_DISABLE_DECODE_CRC32_CHECKS (0)
#endif

// Using unaligned loads and stores causes errors when using UBSan. Jam it off.
#if defined(__has_feature)
	#if __has_feature(undefined_behavior_sanitizer)
		#undef FPNG_USE_UNALIGNED_LOADS
		#define FPNG_USE_UNALIGNED_LOADS (0)
	#endif
#endif

// Set to 0 if your platform doesn't support unaligned 32-bit/64-bit reads/writes. 
#ifndef FPNG_USE_UNALIGNED_LOADS
	#if FPNG_X86_OR_X64_CPU
		// On x86/x64 we default to enabled, for a noticeable perf gain.
		#define FPNG_USE_UNALIGNED_LOADS (1)
	#else
		#define FPNG_USE_UNALIGNED_LOADS (0)
	#endif
#endif

#if defined(_MSC_VER) || defined(__MINGW32__) || FPNG_X86_OR_X64_CPU
	#ifndef __LITTLE_ENDIAN
	#define __LITTLE_ENDIAN 1234
	#endif
	#ifndef __BIG_ENDIAN
	#define __BIG_ENDIAN 4321
	#endif

	// Assume little endian on Windows/x86/x64.
	#define __BYTE_ORDER __LITTLE_ENDIAN
#elif defined(__APPLE__)
	#define __BYTE_ORDER __BYTE_ORDER__
	#define __LITTLE_ENDIAN __LITTLE_ENDIAN__
	#define __BIG_ENDIAN __BIG_ENDIAN__
#else
	// for __BYTE_ORDER (__LITTLE_ENDIAN or __BIG_ENDIAN)
	#include <sys/param.h>

	#ifndef __LITTLE_ENDIAN
	#define __LITTLE_ENDIAN 1234
	#endif
	#ifndef __BIG_ENDIAN
	#define __BIG_ENDIAN 4321
	#endif
#endif

#if !defined(__BYTE_ORDER)
	#error __BYTE_ORDER undefined. Compile with -D__BYTE_ORDER=1234 for little endian or -D__BYTE_ORDER=4321 for big endian.
#endif

namespace fpng
{
	static const int FPNG_FALSE = 0;
	static const uint8_t FPNG_FDEC_VERSION = 0;
	static const uint32_t FPNG_MAX_SUPPORTED_DIM = 1 << 24;

	template <typename S> static inline S maximum(S a, S b) { return (a > b) ? a : b; }
	template <typename S> static inline S minimum(S a, S b) { return (a < b) ? a : b; }

	static inline uint32_t simple_swap32(uint32_t x) { return (x >> 24) | ((x >> 8) & 0x0000FF00) | ((x << 8) & 0x00FF0000) | (x << 24); }
	static inline uint64_t simple_swap64(uint64_t x) { return (((uint64_t)simple_swap32((uint32_t)x)) << 32U) | simple_swap32((uint32_t)(x >> 32U)); }

	static inline uint32_t swap32(uint32_t x)
	{
#if defined(__GNUC__) || defined(__clang__)
		return __builtin_bswap32(x);
#else
		return simple_swap32(x);
#endif
	}

	static inline uint64_t swap64(uint64_t x)
	{
#if defined(__GNUC__) || defined(__clang__)
		return __builtin_bswap64(x);
#else
		return simple_swap64(x);
#endif
	}

#if FPNG_USE_UNALIGNED_LOADS
	#if __BYTE_ORDER == __BIG_ENDIAN
		#define READ_LE32(p) swap32(*reinterpret_cast<const uint32_t *>(p))
		#define WRITE_LE32(p, v) *reinterpret_cast<uint32_t *>(p) = swap32((uint32_t)(v))
		#define WRITE_LE64(p, v) *reinterpret_cast<uint64_t *>(p) = swap64((uint64_t)(v))

		#define READ_BE32(p) *reinterpret_cast<const uint32_t *>(p)
	#else
		#define READ_LE32(p) (*reinterpret_cast<const uint32_t *>(p))
		#define WRITE_LE32(p, v) *reinterpret_cast<uint32_t *>(p) = (uint32_t)(v)
		#define WRITE_LE64(p, v) *reinterpret_cast<uint64_t *>(p) = (uint64_t)(v)

		#define READ_BE32(p) swap32(*reinterpret_cast<const uint32_t *>(p))
	#endif
#else
	// A good compiler should be able to optimize these routines - hopefully. They are crucial for performance.
	static inline uint32_t READ_LE32(const void* p)
	{
		const uint8_t* pBytes = (const uint8_t*)p;
		return ((uint32_t)pBytes[0]) | (((uint32_t)pBytes[1]) << 8U) | (((uint32_t)pBytes[2]) << 16U) | (((uint32_t)pBytes[3]) << 24U);
	}

	static inline uint32_t READ_BE32(const void* p)
	{
		const uint8_t* pBytes = (const uint8_t*)p;
		return ((uint32_t)pBytes[3]) | (((uint32_t)pBytes[2]) << 8U) | (((uint32_t)pBytes[1]) << 16U) | (((uint32_t)pBytes[0]) << 24U);
	}

	static inline void WRITE_LE32(const void* p, uint32_t v)
	{
		uint8_t* pBytes = (uint8_t*)p;
		pBytes[0] = (uint8_t)(v);
		pBytes[1] = (uint8_t)(v >> 8);
		pBytes[2] = (uint8_t)(v >> 16);
		pBytes[3] = (uint8_t)(v >> 24);
	}

	static inline void WRITE_LE64(const void* p, uint64_t v)
	{
		uint8_t* pBytes = (uint8_t*)p;
		pBytes[0] = (uint8_t)(v);
		pBytes[1] = (uint8_t)(v >> 8);
		pBytes[2] = (uint8_t)(v >> 16);
		pBytes[3] = (uint8_t)(v >> 24);
		pBytes[4] = (uint8_t)(v >> 32);
		pBytes[5] = (uint8_t)(v >> 40);
		pBytes[6] = (uint8_t)(v >> 48);
		pBytes[7] = (uint8_t)(v >> 56);
	}
#endif

	// Customized the very common case of reading a 24bpp pixel from memory
	static inline uint32_t READ_RGB_PIXEL(const void* p)
	{
#if FPNG_USE_UNALIGNED_LOADS 
		return READ_LE32(p) & 0xFFFFFF;
#else
		const uint8_t* pBytes = (const uint8_t*)p;
		return ((uint32_t)pBytes[0]) | (((uint32_t)pBytes[1]) << 8U) | (((uint32_t)pBytes[2]) << 16U);
#endif
	}

	// See "Slicing by 4" CRC-32 algorithm here: 
	// https://create.stephan-brumme.com/crc32/

	// Precomputed 4KB of CRC-32 tables
	static const uint32_t g_crc32_4[4][256] = {
	{00, 016701630226, 035603460454, 023102250672, 0733342031, 016032572217, 035130722465, 023631112643, 01666704062, 017167134244, 034065364436, 022764554610, 01155446053, 017654276275, 034756026407, 022057616621, 03555610144, 015254020362, 036356270510, 020457440736, 03266552175, 015567362353, 036465132521, 020364702707, 02333114126, 014432724300, 037530574572, 021231344754, 02400256117, 014301466331, 037203636543, 021502006765,
	07333420310, 011432210136, 032530040744, 024231670562, 07400762321, 011301152107, 032203302775, 024502532553, 06555324372, 010254514154, 033356744726, 025457174500, 06266066343, 010567656165, 033465406717, 025364236531, 04666230254, 012167400072, 031065650600, 027764060426, 04155172265, 012654742043, 031756512631, 027057322417, 05000534236, 013701304010, 030603154662, 026102764444, 05733676207, 013032046021, 030130216653, 026631426475,
	016667040620, 0166670406, 023064420274, 035765210052, 016154302611, 0655532437, 023757762245, 035056152063, 017001744642, 01700174464, 022602324216, 034103514030, 017732406673, 01033236455, 022131066227, 034630656001, 015332650764, 03433060542, 020531230330, 036230400116, 015401512755, 03300322573, 020202172301, 036503742127, 014554154706, 02255764520, 021357534352, 037456304174, 014267216737, 02566426511, 021464676363, 037365046145,
	011554460530, 07255250716, 024357000164, 032456630342, 011267722501, 07566112727, 024464342155, 032365572373, 010332364552, 06433554774, 025531704106, 033230134320, 010401026563, 06300616745, 025202446137, 033503276311, 012001270474, 04700440652, 027602610020, 031103020206, 012732132445, 04033702663, 027131552011, 031630362237, 013667574416, 05166344630, 026064114042, 030765724264, 013154636427, 05655006601, 026757256073, 030056466255,
	035556101440, 023257731666, 0355561014, 016454351232, 035265243471, 023564473657, 0466623025, 016367013203, 034330605422, 022431035604, 01533265076, 017232455250, 034403547413, 022302377635, 01200127047, 017501717261, 036003711504, 020702121722, 03600371150, 015101541376, 036730453535, 020031263713, 03133033161, 015632603347, 037665015566, 021164625740, 02066475132, 014767245314, 037156357557, 021657567771, 02755737103, 014054107325,
	032665521750, 024164311576, 07066141304, 011767771122, 032156663761, 024657053547, 07755203335, 011054433113, 033003225732, 025702415514, 06600645366, 010101075140, 033730167703, 025031757525, 06133507357, 010632337171, 031330331614, 027431501432, 04533751240, 012232161066, 031403073625, 027302643403, 04200413271, 012501223057, 030556435676, 026257205450, 05355055222, 013454665004, 030265777647, 026564147461, 05466317213, 013367527035,
	023331141260, 035430771046, 016532521634, 0233311412, 023402203251, 035303433077, 016201663605, 0500053423, 022557645202, 034256075024, 017354225656, 01455415470, 022264507233, 034565337015, 017467167667, 01366757441, 020664751324, 036165161102, 015067331770, 03766501556, 020157413315, 036656223133, 015754073741, 03055643567, 021002055346, 037703665160, 014601435712, 02100205534, 021731317377, 037030527151, 014132777723, 02633147505,
	024002561170, 032703351356, 011601101524, 07100731702, 024731623141, 032030013367, 011132243515, 07633473733, 025664265112, 033165455334, 010067605546, 06766035760, 025157127123, 033656717305, 010754547577, 06055377751, 027557371034, 031256541212, 012354711460, 04455121646, 027264033005, 031565603223, 012467453451, 04366263677, 026331475056, 030430245270, 013532015402, 05233625624, 026402737067, 030303107241, 013201357433, 05500567615,
	}, { 00,03106630501,06215461202,05313251703,014433142404,017535772105,012626523606,011720313307,031066305010,032160535511,037273764212,034375154713,025455247414,026553477115,023640626616,020746016317,011260411121,012366221420,017075070323,014173640622,05653553525,06755363024,03446132727,0540702226,020206714131,023300124430,026013375333,025115545632,034635656535,037733066034,032420237737,031526407236,
	022541022242,021447612743,024754443040,027652273541,036172160646,035074750347,030367501444,033261331145,013527327252,010421517753,015732746050,016634176551,07114265656,04012455357,01301604454,02207034155,033721433363,030627203662,035534052161,036432662460,027312571767,024214341266,021107110565,022001720064,02747736373,01641106672,04552357171,07454567470,016374674777,015272044276,010161215575,013067425074,
	036036247405,035130477104,030223626607,033325016306,022405305001,021503535500,024610764203,027716154702,07050142415,04156772114,01245523617,02343313316,013463000011,010565630510,015676461213,016770251712,027256656524,024350066025,021043237726,022145407227,033665714120,030763124421,035470375322,036576545623,016230553534,015336363035,010025132736,013123702237,02603411130,01705221431,04416070332,07510640633,
	014577265647,017471455346,012762604445,011664034144,0144327243,03042517742,06351746041,05257176540,025511160657,026417750356,023704501455,020602331154,031122022253,032024612752,037337443051,034231273550,05717674766,06611044267,03502215564,0404425065,011324736362,012222106663,017131357160,014037567461,034771571776,037677341277,032564110574,031462720075,020342433372,023244203673,026157052170,025051662471,
	07340714113,04246124412,01155375311,02053545610,013773656517,010675066016,015566237715,016460407214,036326411103,035220221402,030133070301,033035640600,022715553507,021613363006,024500132705,027406702204,016120305032,015026535533,010335764230,013233154731,02513247436,01415477137,04706626634,07600016335,027146000022,024040630523,021353461220,022255251721,033575142426,030473772127,035760523624,036666313325,
	025601736351,026707106650,023414357153,020512567452,031232674755,032334044254,037027215557,034121425056,014667433341,017761203640,012472052143,011574662442,0254571745,03352341244,06041110547,05147720046,034461327270,037567517771,032674746072,031772176573,020052265674,023154455375,026247604476,025341034177,05407022260,06501612761,03612443062,0714273563,011034160664,012132750365,017221501466,014327331167,
	031376553516,032270363017,037163132714,034065702215,025745411112,026643221413,023550070310,020456640611,0310656506,03216066007,06105237704,05003407205,014723714102,017625124403,012536375300,011430545601,020116142437,023010772136,026303523635,025205313334,034525000033,037423630532,032730461231,031636251730,011170247427,012076477126,017365626625,014263016324,05543305023,06445535522,03756764221,0650154720,
	013637571754,010731341255,015422110556,016524720057,07204433350,04302203651,01011052152,02117662453,022651674744,021757044245,024444215546,027542425047,036262736340,035364106641,030077357142,033171567443,02457160675,01551750374,04642501477,07744331176,016064022271,015162612770,010271443073,013377273572,033431265665,030537455364,035624604467,036722034166,027002327261,024104517760,021217746063,022311176562,
	}, { 00,0160465067,0341152156,0221537131,0702324334,0662741353,0443276262,0523613205,01604650670,01764235617,01545702726,01425367741,01106574544,01066111523,01247426412,01327043475,03411521560,03571144507,03750473436,03630016451,03313605654,03273260633,03052757702,03132332765,02215371310,02375714377,02154223246,02034646221,02517055024,02477430043,02656107172,02736562115,
	07023243340,07143626327,07362311216,07202774271,07721167074,07641502013,07460035122,07500450145,06627413530,06747076557,06566541466,06406124401,06125737604,06045352663,06264665752,06304200735,04432762620,04552307647,04773630776,04613255711,04330446514,04250023573,04071514442,04111171425,05236132050,05356557037,05177060106,05017405161,05534216364,05454673303,05675344232,05715721255,
	016046506700,016126163767,016307454656,016267031631,016744622434,016624247453,016405770562,016565315505,017642356170,017722733117,017503204026,017463661041,017140072244,017020417223,017201120312,017361545375,015457027260,015537442207,015716175336,015676510351,015355303154,015235766133,015014251002,015174634065,014253677410,014333212477,014112725546,014072340521,014551553724,014431136743,014610401672,014770064615,
	011065745440,011105320427,011324617516,011244272571,011767461774,011607004713,011426533622,011546156645,010661115230,010701570257,010520047366,010440422301,010163231104,010003654163,010222363052,010342706035,012474264120,012514601147,012735336076,012655753011,012376140214,012216525273,012037012342,012157477325,013270434750,013310051737,013131566606,013051103661,013572710464,013412375403,013633642532,013753227555,
	034115215600,034075670667,034254347756,034334722731,034617131534,034777554553,034556063462,034436406405,035711445070,035671020017,035450517126,035530172141,035013761344,035173304323,035352633212,035232256275,037504734360,037464351307,037645666236,037725203251,037206410054,037366075033,037147542102,037027127165,036300164510,036260501577,036041036446,036121453421,036402240624,036562625643,036743312772,036623777715,
	033136056540,033056433527,033277104416,033317561471,033634372674,033754717613,033575220722,033415645745,032732606330,032652263357,032473754266,032513331201,032030522004,032150147063,032371470152,032211015135,030527577020,030447112047,030666425176,030706040111,030225653314,030345236373,030164701242,030004364225,031323327650,031243742637,031062275706,031102610761,031421003564,031541466503,031760151432,031600534455,
	022153713100,022033376167,022212641056,022372224031,022651437234,022731052253,022510565362,022470100305,023757143770,023637526717,023416011626,023576474641,023055267444,023135602423,023314335512,023274750575,021542232460,021422657407,021603360536,021763705551,021240116754,021320573733,021101044602,021061421665,020346462210,020226007277,020007530346,020167155321,020444746124,020524323143,020705614072,020665271015,
	025170550240,025010135227,025231402316,025351067371,025672674174,025712211113,025533726022,025453343045,024774300430,024614765457,024435252566,024555637501,024076024704,024116441763,024337176652,024257513635,026561071720,026401414747,026620123676,026740546611,026263355414,026303730473,026122207542,026042662525,027365621150,027205244137,027024773006,027144316061,027467505264,027507160203,027726457332,027646032355,
	}, { 00,027057063545,025202344213,02255327756,021730513527,06767570062,04532657734,023565634271,030555024357,017502047612,015757360144,032700303401,011265537670,036232554335,034067673463,013030610126,012006253637,035051230372,037204117424,010253174161,033736740310,014761723655,016534404103,031563467446,022553277560,05504214025,07751133773,020706150236,03263764047,024234707502,026061420254,01036443711,
	024014527476,03043544133,01216663665,026241600320,05724034151,022773057414,020526370342,07571313607,014541503721,033516560264,031743647532,016714624077,035271010206,012226073743,010073354015,037024337550,036012774241,011045717704,013210430052,034247453517,017722267766,030775204223,032520123575,015577140030,06547750116,021510733453,023745414305,04712477640,027277243431,0220220174,02075107622,025022164367,
	023305054075,04352037530,06107310266,021150373723,02435547552,025462524017,027637603741,0660660204,013650070322,034607013667,036452334131,011405357474,032160563605,015137500340,017362627416,030335644153,031303207642,016354264307,014101143451,033156120114,010433714365,037464777620,035631450176,012666433433,01656223515,026601240050,024454167706,03403104243,020166730032,07131753577,05364474221,022333417764,
	07311573403,020346510146,022113637610,05144654355,026421060124,01476003461,03623324337,024674347672,037644557754,010613534211,012446613547,035411670002,016174044273,031123027736,033376300060,014321363525,015317720234,032340743771,030115464027,017142407562,034427233713,013470250256,011625177500,036672114045,025642704163,02615767426,0440440370,027417423635,04172217444,023125274101,021370153657,06327130312,
	035526333073,012571350536,010724077260,037773014725,014216620554,033241643011,031014564747,016043507202,05073317324,022024374661,020271053137,07226030472,024743604603,03714667346,01541540410,026516523155,027520160644,0577103301,02722224457,025775247112,06210473363,021247410626,023012737170,04045754435,017075144513,030022127056,032277200700,015220263245,036745457034,011712434571,013547713227,034510770762,
	011532614405,036565677140,034730550616,013767533353,030202307122,017255364467,015000043331,032057020674,021067630752,06030653217,04265574541,023232517004,0757323275,027700340730,025555067066,02502004523,03534447232,024563424777,026736703021,01761760564,022204154715,05253137250,07006210506,020051273043,033061463165,014036400420,016263727376,031234744633,012751170442,035706113107,037553234651,010504257314,
	016623367006,031674304543,033421023215,014476040750,037113674521,010144617064,012311530732,035346553277,026376343351,01321320614,03174007142,024123064407,07446650676,020411633333,022644514465,05613577120,04625134631,023672157374,021427270422,06470213167,025115427316,02142444653,0317763105,027340700440,034370110566,013327173023,011172254775,036125237230,015440403041,032417460504,030642747252,017615724717,
	032637640470,015660623135,017435504663,030462567326,013107353157,034150330412,036305017344,011352074601,02362664727,025335607262,027160520534,0137543071,023452377200,04405314745,06650033013,021607050556,020631413247,07666470702,05433757054,022464734511,01101100760,026156163225,024303244573,03354227036,010364437110,037333454455,035166773303,012131710646,031454124437,016403147172,014656260624,033601203361,
	} };

	static uint32_t crc32_slice_by_4(const void* pData, size_t data_len, uint32_t cur_crc32 = 0)
	{
		uint32_t crc = ~cur_crc32;
		const uint32_t* pData32 = static_cast<const uint32_t*>(pData);

		for (; data_len >= sizeof(uint32_t); ++pData32, data_len -= 4)
		{
			uint32_t v = READ_LE32(pData32) ^ crc;
			crc = g_crc32_4[0][v >> 24] ^ g_crc32_4[1][(v >> 16) & 0xFF] ^ g_crc32_4[2][(v >> 8) & 0xFF] ^ g_crc32_4[3][v & 0xFF];
		}

		for (const uint8_t* pData8 = reinterpret_cast<const uint8_t*>(pData32); data_len; --data_len)
			crc = (crc >> 8) ^ g_crc32_4[0][(crc & 0xFF) ^ *pData8++];

		return ~crc;
	}

#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE 
	// See Fast CRC Computation for Generic Polynomials Using PCLMULQDQ Instruction":
	// https://www.intel.com/content/dam/www/public/us/en/documents/white-papers/fast-crc-computation-generic-polynomials-pclmulqdq-paper.pdf
	// Requires PCLMUL and SSE 4.1. This function skips Step 1 (fold by 4) for simplicity/less code.
	static uint32_t crc32_pclmul(const uint8_t* p, size_t size, uint32_t crc)
	{
		assert(size >= 16);

		// See page 22 (bit reflected constants for gzip)
#ifdef _MSC_VER
		static const uint64_t __declspec(align(16)) 
#else
		static const uint64_t __attribute__((aligned(16)))
#endif
			s_u[2] = { 0x1DB710641, 0x1F7011641 }, s_k5k0[2] = { 0x163CD6124, 0 }, s_k3k4[2] = { 0x1751997D0, 0xCCAA009E };

		// Load first 16 bytes, apply initial CRC32
		__m128i b = _mm_xor_si128(_mm_cvtsi32_si128(~crc), _mm_loadu_si128(reinterpret_cast<const __m128i*>(p)));

		// We're skipping directly to Step 2 page 12 - iteratively folding by 1 (by 4 is overkill for our needs)
		const __m128i k3k4 = _mm_load_si128(reinterpret_cast<const __m128i*>(s_k3k4));

		for (size -= 16, p += 16; size >= 16; size -= 16, p += 16)
			b = _mm_xor_si128(_mm_xor_si128(_mm_clmulepi64_si128(b, k3k4, 17), _mm_loadu_si128(reinterpret_cast<const __m128i*>(p))), _mm_clmulepi64_si128(b, k3k4, 0));

		// Final stages: fold to 64-bits, 32-bit Barrett reduction
		const __m128i z = _mm_set_epi32(0, ~0, 0, ~0), u = _mm_load_si128(reinterpret_cast<const __m128i*>(s_u));
		b = _mm_xor_si128(_mm_srli_si128(b, 8), _mm_clmulepi64_si128(b, k3k4, 16));
		b = _mm_xor_si128(_mm_clmulepi64_si128(_mm_and_si128(b, z), _mm_loadl_epi64(reinterpret_cast<const __m128i*>(s_k5k0)), 0), _mm_srli_si128(b, 4));
		return ~_mm_extract_epi32(_mm_xor_si128(b, _mm_clmulepi64_si128(_mm_and_si128(_mm_clmulepi64_si128(_mm_and_si128(b, z), u, 16), z), u, 0)), 1);
	}

	static uint32_t crc32_sse41_simd(const unsigned char* buf, size_t len, uint32_t prev_crc32)
	{
		if (len < 16)
			return crc32_slice_by_4(buf, len, prev_crc32);

		uint32_t simd_len = len & ~15;
		uint32_t c = crc32_pclmul(buf, simd_len, prev_crc32);
		return crc32_slice_by_4(buf + simd_len, len - simd_len, c);
	}
#endif

#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE 

#ifndef _MSC_VER
	static void do_cpuid(uint32_t eax, uint32_t ecx, uint32_t* regs)
	{
		uint32_t ebx = 0, edx = 0;

#if defined(__PIC__) && defined(__i386__)
		__asm__("movl %%ebx, %%edi;"
			"cpuid;"
			"xchgl %%ebx, %%edi;"
			: "=D"(ebx), "+a"(eax), "+c"(ecx), "=d"(edx));
#else
		__asm__("cpuid;" : "+b"(ebx), "+a"(eax), "+c"(ecx), "=d"(edx));
#endif

		regs[0] = eax; regs[1] = ebx; regs[2] = ecx; regs[3] = edx;
	}
#endif

	struct cpu_info
	{
		cpu_info() { memset(this, 0, sizeof(*this)); }

		bool m_initialized, m_has_fpu, m_has_mmx, m_has_sse, m_has_sse2, m_has_sse3, m_has_ssse3, m_has_sse41, m_has_sse42, m_has_avx, m_has_avx2, m_has_pclmulqdq;
				
		void init()
		{
			if (m_initialized)
				return;

			int regs[4];

#ifdef _MSC_VER
			__cpuid(regs, 0);
#else
			do_cpuid(0, 0, (uint32_t*)regs);
#endif

			const uint32_t max_eax = regs[0];
			if (max_eax >= 1U)
			{
#ifdef _MSC_VER
				__cpuid(regs, 1);
#else
				do_cpuid(1, 0, (uint32_t*)regs);
#endif
				extract_x86_flags(regs[2], regs[3]);
			}

			if (max_eax >= 7U)
			{
#ifdef _MSC_VER
				__cpuidex(regs, 7, 0);
#else
				do_cpuid(7, 0, (uint32_t*)regs);
#endif
				extract_x86_extended_flags(regs[1]);
			}

			m_initialized = true;
		}

		bool can_use_sse41() const { return m_has_sse && m_has_sse2 && m_has_sse3 && m_has_ssse3 && m_has_sse41; }
		bool can_use_pclmul() const	{ return m_has_pclmulqdq && can_use_sse41(); }

	private:
		void extract_x86_flags(uint32_t ecx, uint32_t edx)
		{
			m_has_fpu = (edx & (1 << 0)) != 0;	m_has_mmx = (edx & (1 << 23)) != 0;	m_has_sse = (edx & (1 << 25)) != 0; m_has_sse2 = (edx & (1 << 26)) != 0;
			m_has_sse3 = (ecx & (1 << 0)) != 0; m_has_ssse3 = (ecx & (1 << 9)) != 0; m_has_sse41 = (ecx & (1 << 19)) != 0; m_has_sse42 = (ecx & (1 << 20)) != 0;
			m_has_pclmulqdq = (ecx & (1 << 1)) != 0; m_has_avx = (ecx & (1 << 28)) != 0;
		}

		void extract_x86_extended_flags(uint32_t ebx) { m_has_avx2 = (ebx & (1 << 5)) != 0; }
	};

	cpu_info g_cpu_info;
		
	void fpng_init()
	{
		g_cpu_info.init();
	}
#else
	void fpng_init()
	{
	}
#endif

	bool fpng_cpu_supports_sse41()
	{
#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE 
		assert(g_cpu_info.m_initialized);
		return g_cpu_info.can_use_sse41();
#else
		return false;
#endif
	}

	uint32_t fpng_crc32(const void* pData, size_t size, uint32_t prev_crc32)
	{
#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE 
		if (g_cpu_info.can_use_pclmul())
			return crc32_sse41_simd(static_cast<const uint8_t *>(pData), size, prev_crc32);
#endif

		return crc32_slice_by_4(pData, size, prev_crc32);
	}

#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE 
	// See "Fast Computation of Adler32 Checksums":
	// https://www.intel.com/content/www/us/en/developer/articles/technical/fast-computation-of-adler32-checksums.html
	// SSE 4.1, 16 bytes per iteration
	static uint32_t adler32_sse_16(const uint8_t* p, size_t len, uint32_t initial)
	{
		uint32_t s1 = initial & 0xFFFF, s2 = initial >> 16;
		const uint32_t K = 65521;

		while (len >= 16)
		{
			__m128i a = _mm_setr_epi32(s1, 0, 0, 0), b = _mm_setzero_si128(), c = _mm_setzero_si128(), d = _mm_setzero_si128(), 
				e = _mm_setzero_si128(), f = _mm_setzero_si128(), g = _mm_setzero_si128(), h = _mm_setzero_si128();

			const size_t n = minimum<size_t>(len >> 4, 5552);

			for (size_t i = 0; i < n; i++)
			{
				const __m128i v = _mm_loadu_si128((const __m128i*)(p + i * 16));
				a = _mm_add_epi32(a, _mm_cvtepu8_epi32(_mm_shuffle_epi32(v, _MM_SHUFFLE(0, 0, 0, 0)))); b = _mm_add_epi32(b, a);
				c = _mm_add_epi32(c, _mm_cvtepu8_epi32(_mm_shuffle_epi32(v, _MM_SHUFFLE(1, 1, 1, 1)))); d = _mm_add_epi32(d, c);
				e = _mm_add_epi32(e, _mm_cvtepu8_epi32(_mm_shuffle_epi32(v, _MM_SHUFFLE(2, 2, 2, 2)))); f = _mm_add_epi32(f, e);
				g = _mm_add_epi32(g, _mm_cvtepu8_epi32(_mm_shuffle_epi32(v, _MM_SHUFFLE(3, 3, 3, 3)))); h = _mm_add_epi32(h, g);
			}

			uint32_t sa[16], sb[16];
			_mm_storeu_si128((__m128i*)sa, a); _mm_storeu_si128((__m128i*)(sa + 4), c);
			_mm_storeu_si128((__m128i*)sb, b); _mm_storeu_si128((__m128i*)(sb + 4), d);
			_mm_storeu_si128((__m128i*)(sa + 8), e); _mm_storeu_si128((__m128i*)(sa + 12), g);
			_mm_storeu_si128((__m128i*)(sb + 8), f); _mm_storeu_si128((__m128i*)(sb + 12), h);

			// This could be vectorized, but it's only executed every 5552*16 iterations.
			uint64_t vs1 = 0;
			for (uint32_t i = 0; i < 16; i++)
				vs1 += sa[i];

			uint64_t vs2_a = 0;
			for (uint32_t i = 0; i < 16; i++)
				vs2_a += sa[i] * (uint64_t)i;
			uint64_t vs2_b = 0;
			for (uint32_t i = 0; i < 16; i++)
				vs2_b += sb[i];
			vs2_b *= 16U;
			uint64_t vs2 = vs2_b - vs2_a + s2;

			s1 = (uint32_t)(vs1 % K);
			s2 = (uint32_t)(vs2 % K);

			p += n * 16;
			len -= n * 16;
		}

		for (; len; len--)
		{
			s1 += *p++;
			s2 += s1;
		}

		return (s1 % K) | ((s2 % K) << 16);
	}
#endif

	static uint32_t fpng_adler32_scalar(const uint8_t* ptr, size_t buf_len, uint32_t adler)
	{
		uint32_t i, s1 = (uint32_t)(adler & 0xffff), s2 = (uint32_t)(adler >> 16); uint32_t block_len = (uint32_t)(buf_len % 5552);
		if (!ptr) return FPNG_ADLER32_INIT;
		while (buf_len) {
			for (i = 0; i + 7 < block_len; i += 8, ptr += 8) {
				s1 += ptr[0], s2 += s1; s1 += ptr[1], s2 += s1; s1 += ptr[2], s2 += s1; s1 += ptr[3], s2 += s1;
				s1 += ptr[4], s2 += s1; s1 += ptr[5], s2 += s1; s1 += ptr[6], s2 += s1; s1 += ptr[7], s2 += s1;
			}
			for (; i < block_len; ++i) s1 += *ptr++, s2 += s1;
			s1 %= 65521U, s2 %= 65521U; buf_len -= block_len; block_len = 5552;
		}
		return (s2 << 16) + s1;
	}

	uint32_t fpng_adler32(const void* pData, size_t size, uint32_t adler)
	{
#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE 
		if (g_cpu_info.can_use_sse41())
			return adler32_sse_16((const uint8_t*)pData, size, adler);
#endif
		return fpng_adler32_scalar((const uint8_t*)pData, size, adler);
	}

	// Ensure we've been configured for endianness correctly.
	static inline bool endian_check()
	{
		uint32_t endian_check = 0;
		WRITE_LE32(&endian_check, 0x1234ABCD);
		const uint32_t first_byte = reinterpret_cast<const uint8_t*>(&endian_check)[0];
		return first_byte == 0xCD;
	}
		
	static const uint16_t g_defl_len_sym[256] = {
	  257,258,259,260,261,262,263,264,265,265,266,266,267,267,268,268,269,269,269,269,270,270,270,270,271,271,271,271,272,272,272,272,
	  273,273,273,273,273,273,273,273,274,274,274,274,274,274,274,274,275,275,275,275,275,275,275,275,276,276,276,276,276,276,276,276,
	  277,277,277,277,277,277,277,277,277,277,277,277,277,277,277,277,278,278,278,278,278,278,278,278,278,278,278,278,278,278,278,278,
	  279,279,279,279,279,279,279,279,279,279,279,279,279,279,279,279,280,280,280,280,280,280,280,280,280,280,280,280,280,280,280,280,
	  281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,281,
	  282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,282,
	  283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,283,
	  284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,284,285 };

	static const uint8_t g_defl_len_extra[256] = {
	  0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
	  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
	  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0 };

	static const uint8_t g_defl_small_dist_sym[512] = {
	  0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,
	  11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,
	  13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
	  14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,
	  14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	  15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,
	  16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
	  16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
	  16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,
	  17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,
	  17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,
	  17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17 };
		
	static const uint32_t g_bitmasks[17] = { 0x0000, 0x0001, 0x0003, 0x0007, 0x000F, 0x001F, 0x003F, 0x007F, 0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF };

	// Huffman tables generated by fpng_test -t @filelist.txt. Total alpha files : 1440, Total opaque files : 5627.
	// Feel free to retrain the encoder on your opaque/alpha PNG files by setting FPNG_TRAIN_HUFFMAN_TABLES and running fpng_test with the -t option.
	static const uint8_t g_dyn_huff_3[] = {
	120, 1, 237, 195, 3, 176, 110, 89, 122, 128, 225, 247, 251, 214, 218, 248, 113, 124, 173, 190, 109, 12, 50, 201, 196, 182, 109, 219, 182, 109, 219, 182,
	109, 219, 201, 36, 147, 153, 105, 235, 246, 53, 142, 207, 143, 141, 181, 214, 151, 93, 117, 170, 78, 117, 117, 58, 206, 77, 210, 217, 169, 122 };
	const uint32_t DYN_HUFF_3_BITBUF = 30, DYN_HUFF_3_BITBUF_SIZE = 7;
	static const struct { uint8_t m_code_size; uint16_t m_code; } g_dyn_huff_3_codes[288] = {
	{2,0},{4,2},{4,10},{5,14},{5,30},{6,25},{6,57},{6,5},{6,37},{7,3},{7,67},{7,35},{7,99},{8,11},{8,139},{8,75},{8,203},{8,43},{8,171},{8,107},{9,135},{9,391},{9,71},{9,327},{9,199},{9,455},{9,39},{9,295},{9,167},{9,423},{9,103},{10,183},
	{9,359},{10,695},{10,439},{10,951},{10,119},{10,631},{10,375},{10,887},{10,247},{10,759},{10,503},{11,975},{11,1999},{11,47},{11,1071},{12,1199},{11,559},{12,3247},{12,687},{11,1583},{12,2735},{12,1711},{12,3759},{12,431},{12,2479},{12,1455},{12,3503},{12,943},{12,2991},{12,1967},{12,4015},{12,111},
	{12,2159},{12,1135},{12,3183},{12,623},{12,2671},{12,1647},{12,3695},{12,367},{12,2415},{12,1391},{12,3439},{12,879},{12,2927},{12,1903},{12,3951},{12,239},{12,2287},{12,1263},{12,3311},{12,751},{12,2799},{12,1775},{12,3823},{12,495},{12,2543},{12,1519},{12,3567},{12,1007},{12,3055},{12,2031},{12,4079},{12,31},
	{12,2079},{12,1055},{12,3103},{12,543},{12,2591},{12,1567},{12,3615},{12,287},{12,2335},{12,1311},{12,3359},{12,799},{12,2847},{12,1823},{12,3871},{12,159},{12,2207},{12,1183},{12,3231},{12,671},{12,2719},{12,1695},{12,3743},{12,415},{12,2463},{12,1439},{12,3487},{12,927},{12,2975},{12,1951},{12,3999},{12,95},
	{12,2143},{12,1119},{12,3167},{12,607},{12,2655},{12,1631},{12,3679},{12,351},{12,2399},{12,1375},{12,3423},{12,863},{12,2911},{12,1887},{12,3935},{12,223},{12,2271},{12,1247},{12,3295},{12,735},{12,2783},{12,1759},{12,3807},{12,479},{12,2527},{12,1503},{12,3551},{12,991},{12,3039},{12,2015},{12,4063},{12,63},
	{12,2111},{12,1087},{12,3135},{12,575},{12,2623},{12,1599},{12,3647},{12,319},{12,2367},{12,1343},{12,3391},{12,831},{12,2879},{12,1855},{12,3903},{12,191},{12,2239},{12,1215},{12,3263},{12,703},{12,2751},{12,1727},{12,3775},{12,447},{12,2495},{12,1471},{12,3519},{12,959},{12,3007},{12,1983},{12,4031},{12,127},
	{12,2175},{12,1151},{12,3199},{12,639},{12,2687},{12,1663},{12,3711},{12,383},{12,2431},{12,1407},{12,3455},{12,895},{12,2943},{11,303},{12,1919},{12,3967},{11,1327},{12,255},{11,815},{11,1839},{11,175},{10,1015},{10,15},{10,527},{10,271},{10,783},{10,143},{10,655},{10,399},{10,911},{10,79},{10,591},
	{9,231},{10,335},{9,487},{9,23},{9,279},{9,151},{9,407},{9,87},{9,343},{9,215},{9,471},{9,55},{8,235},{8,27},{8,155},{8,91},{8,219},{8,59},{8,187},{8,123},{7,19},{7,83},{7,51},{7,115},{6,21},{6,53},{6,13},{6,45},{5,1},{5,17},{5,9},{4,6},
	{12,2303},{6,29},{0,0},{0,0},{8,251},{0,0},{0,0},{8,7},{0,0},{10,847},{0,0},{10,207},{12,1279},{10,719},{12,3327},{12,767},{12,2815},{12,1791},{12,3839},{12,511},{12,2559},{12,1535},{9,311},{12,3583},{12,1023},{12,3071},{10,463},{12,2047},{6,61},{12,4095},{0,0},{0,0}
	};

	static const uint8_t g_dyn_huff_4[] = {
	120, 1, 229, 196, 99, 180, 37, 103, 218, 128, 225, 251, 121, 171, 106, 243, 216, 231, 180, 109, 196, 182, 51, 51, 73, 6, 201, 216, 182, 109, 219, 182,
	17, 140, 98, 219, 102, 219, 60, 125, 172, 205, 170, 122, 159, 111, 213, 143, 179, 214, 94, 189, 58, 153, 104, 166, 103, 190, 247, 199, 117 };
	const uint32_t DYN_HUFF_4_BITBUF = 1, DYN_HUFF_4_BITBUF_SIZE = 2;
	static const struct { uint8_t m_code_size; uint16_t m_code; } g_dyn_huff_4_codes[288] = {
	{2,0},{4,2},{5,6},{6,30},{6,62},{6,1},{7,41},{7,105},{7,25},{7,89},{7,57},{7,121},{8,117},{8,245},{8,13},{8,141},{8,77},{8,205},{8,45},{8,173},{8,109},{8,237},{8,29},{8,157},{8,93},{8,221},{8,61},{9,83},{9,339},{9,211},{9,467},{9,51},
	{9,307},{9,179},{9,435},{9,115},{9,371},{9,243},{9,499},{9,11},{9,267},{9,139},{9,395},{9,75},{9,331},{9,203},{9,459},{9,43},{9,299},{10,7},{10,519},{10,263},{10,775},{10,135},{10,647},{10,391},{10,903},{10,71},{10,583},{10,327},{10,839},{10,199},{10,711},{10,455},
	{10,967},{10,39},{10,551},{10,295},{10,807},{10,167},{10,679},{10,423},{10,935},{10,103},{10,615},{11,463},{11,1487},{11,975},{10,359},{10,871},{10,231},{11,1999},{11,47},{11,1071},{11,559},{10,743},{10,487},{11,1583},{11,303},{11,1327},{11,815},{11,1839},{11,175},{11,1199},{11,687},{11,1711},
	{11,431},{11,1455},{11,943},{11,1967},{11,111},{11,1135},{11,623},{11,1647},{11,367},{11,1391},{11,879},{11,1903},{11,239},{11,1263},{11,751},{11,1775},{11,495},{11,1519},{11,1007},{11,2031},{11,31},{11,1055},{11,543},{11,1567},{11,287},{11,1311},{11,799},{11,1823},{11,159},{11,1183},{11,671},{11,1695},
	{11,415},{11,1439},{11,927},{11,1951},{11,95},{11,1119},{11,607},{11,1631},{11,351},{11,1375},{11,863},{11,1887},{11,223},{11,1247},{11,735},{11,1759},{11,479},{11,1503},{11,991},{11,2015},{11,63},{11,1087},{11,575},{11,1599},{11,319},{11,1343},{11,831},{11,1855},{11,191},{11,1215},{11,703},{11,1727},
	{11,447},{11,1471},{11,959},{11,1983},{11,127},{11,1151},{11,639},{11,1663},{11,383},{10,999},{10,23},{10,535},{10,279},{11,1407},{11,895},{11,1919},{11,255},{11,1279},{10,791},{10,151},{10,663},{10,407},{10,919},{10,87},{10,599},{10,343},{10,855},{10,215},{10,727},{10,471},{10,983},{10,55},
	{10,567},{10,311},{10,823},{10,183},{10,695},{10,439},{10,951},{10,119},{10,631},{10,375},{10,887},{10,247},{10,759},{10,503},{10,1015},{10,15},{10,527},{10,271},{10,783},{10,143},{10,655},{10,399},{9,171},{9,427},{9,107},{9,363},{9,235},{9,491},{9,27},{9,283},{9,155},{9,411},
	{9,91},{9,347},{9,219},{9,475},{9,59},{9,315},{9,187},{9,443},{8,189},{9,123},{8,125},{8,253},{8,3},{8,131},{8,67},{8,195},{8,35},{8,163},{8,99},{8,227},{8,19},{7,5},{7,69},{7,37},{7,101},{7,21},{7,85},{6,33},{6,17},{6,49},{5,22},{4,10},
	{12,2047},{0,0},{6,9},{0,0},{0,0},{0,0},{8,147},{0,0},{0,0},{7,53},{0,0},{9,379},{0,0},{9,251},{10,911},{10,79},{11,767},{10,591},{10,335},{10,847},{10,207},{10,719},{11,1791},{11,511},{9,507},{11,1535},{11,1023},{12,4095},{5,14},{0,0},{0,0},{0,0}
	};

#define PUT_BITS(bb, ll) do { uint32_t b = bb, l = ll; assert((l) >= 0 && (l) <= 16); assert((b) < (1ULL << (l))); bit_buf |= (((uint64_t)(b)) << bit_buf_size); bit_buf_size += (l); assert(bit_buf_size <= 64); } while(0)
#define PUT_BITS_CZ(bb, ll) do { uint32_t b = bb, l = ll; assert((l) >= 1 && (l) <= 16); assert((b) < (1ULL << (l))); bit_buf |= (((uint64_t)(b)) << bit_buf_size); bit_buf_size += (l); assert(bit_buf_size <= 64); } while(0)

#define PUT_BITS_FLUSH do { \
	if ((dst_ofs + 8) > dst_buf_size) \
		return 0; \
	WRITE_LE64(pDst + dst_ofs, bit_buf); \
	uint32_t bits_to_shift = bit_buf_size & ~7; \
	dst_ofs += (bits_to_shift >> 3); \
	assert(bits_to_shift < 64); \
	bit_buf = bit_buf >> bits_to_shift; \
	bit_buf_size -= bits_to_shift; \
} while(0)

#define PUT_BITS_FORCE_FLUSH do { \
	while (bit_buf_size > 0) \
	{ \
		if ((dst_ofs + 1) > dst_buf_size) \
			return 0; \
		*(uint8_t*)(pDst + dst_ofs) = (uint8_t)bit_buf; \
		dst_ofs++; \
		bit_buf >>= 8; \
		bit_buf_size -= 8; \
	} \
} while(0)

	enum
	{
		DEFL_MAX_HUFF_TABLES = 3,
		DEFL_MAX_HUFF_SYMBOLS = 288,	
		DEFL_MAX_HUFF_SYMBOLS_0 = 288,	
		DEFL_MAX_HUFF_SYMBOLS_1 = 32,
		DEFL_MAX_HUFF_SYMBOLS_2 = 19,
		DEFL_LZ_DICT_SIZE = 32768,
		DEFL_LZ_DICT_SIZE_MASK = DEFL_LZ_DICT_SIZE - 1,
		DEFL_MIN_MATCH_LEN = 3,
		DEFL_MAX_MATCH_LEN = 258
	};

#if FPNG_TRAIN_HUFFMAN_TABLES
	uint64_t g_huff_counts[HUFF_COUNTS_SIZE];
#endif

	struct defl_huff
	{
		uint16_t m_huff_count[DEFL_MAX_HUFF_TABLES][DEFL_MAX_HUFF_SYMBOLS];
		uint16_t m_huff_codes[DEFL_MAX_HUFF_TABLES][DEFL_MAX_HUFF_SYMBOLS];
		uint8_t m_huff_code_sizes[DEFL_MAX_HUFF_TABLES][DEFL_MAX_HUFF_SYMBOLS];
	};

	struct defl_sym_freq
	{
		uint16_t m_key;
		uint16_t m_sym_index;
	};

#define DEFL_CLEAR_OBJ(obj) memset(&(obj), 0, sizeof(obj))

	static defl_sym_freq* defl_radix_sort_syms(uint32_t num_syms, defl_sym_freq* pSyms0, defl_sym_freq* pSyms1)
	{
		uint32_t total_passes = 2, pass_shift, pass, i, hist[256 * 2]; defl_sym_freq* pCur_syms = pSyms0, * pNew_syms = pSyms1; DEFL_CLEAR_OBJ(hist);
		for (i = 0; i < num_syms; i++) { uint32_t freq = pSyms0[i].m_key; hist[freq & 0xFF]++; hist[256 + ((freq >> 8) & 0xFF)]++; }
		while ((total_passes > 1) && (num_syms == hist[(total_passes - 1) * 256])) total_passes--;
		for (pass_shift = 0, pass = 0; pass < total_passes; pass++, pass_shift += 8)
		{
			const uint32_t* pHist = &hist[pass << 8];
			uint32_t offsets[256], cur_ofs = 0;
			for (i = 0; i < 256; i++) { offsets[i] = cur_ofs; cur_ofs += pHist[i]; }
			for (i = 0; i < num_syms; i++) pNew_syms[offsets[(pCur_syms[i].m_key >> pass_shift) & 0xFF]++] = pCur_syms[i];
			{ defl_sym_freq* t = pCur_syms; pCur_syms = pNew_syms; pNew_syms = t; }
		}
		return pCur_syms;
	}

	// defl_calculate_minimum_redundancy() originally written by: Alistair Moffat, alistair@cs.mu.oz.au, Jyrki Katajainen, jyrki@diku.dk, November 1996.
	static void defl_calculate_minimum_redundancy(defl_sym_freq* A, int n)
	{
		int root, leaf, next, avbl, used, dpth;
		if (n == 0) return; else if (n == 1) { A[0].m_key = 1; return; }
		A[0].m_key += A[1].m_key; root = 0; leaf = 2;
		for (next = 1; next < n - 1; next++)
		{
			if (leaf >= n || A[root].m_key < A[leaf].m_key) { A[next].m_key = A[root].m_key; A[root++].m_key = (uint16_t)next; }
			else A[next].m_key = A[leaf++].m_key;
			if (leaf >= n || (root < next && A[root].m_key < A[leaf].m_key)) { A[next].m_key = (uint16_t)(A[next].m_key + A[root].m_key); A[root++].m_key = (uint16_t)next; }
			else A[next].m_key = (uint16_t)(A[next].m_key + A[leaf++].m_key);
		}
		A[n - 2].m_key = 0; for (next = n - 3; next >= 0; next--) A[next].m_key = A[A[next].m_key].m_key + 1;
		avbl = 1; used = dpth = 0; root = n - 2; next = n - 1;
		while (avbl > 0)
		{
			while (root >= 0 && (int)A[root].m_key == dpth) { used++; root--; }
			while (avbl > used) { A[next--].m_key = (uint16_t)(dpth); avbl--; }
			avbl = 2 * used; dpth++; used = 0;
		}
	}

	// Limits canonical Huffman code table's max code size.
	enum { DEFL_MAX_SUPPORTED_HUFF_CODESIZE = 32 };
	static void defl_huffman_enforce_max_code_size(int* pNum_codes, int code_list_len, int max_code_size)
	{
		int i; uint32_t total = 0; if (code_list_len <= 1) return;
		for (i = max_code_size + 1; i <= DEFL_MAX_SUPPORTED_HUFF_CODESIZE; i++) pNum_codes[max_code_size] += pNum_codes[i];
		for (i = max_code_size; i > 0; i--) total += (((uint32_t)pNum_codes[i]) << (max_code_size - i));
		while (total != (1UL << max_code_size))
		{
			pNum_codes[max_code_size]--;
			for (i = max_code_size - 1; i > 0; i--) if (pNum_codes[i]) { pNum_codes[i]--; pNum_codes[i + 1] += 2; break; }
			total--;
		}
	}

	static void defl_optimize_huffman_table(defl_huff* d, int table_num, int table_len, int code_size_limit, int static_table)
	{
		int i, j, l, num_codes[1 + DEFL_MAX_SUPPORTED_HUFF_CODESIZE]; uint32_t next_code[DEFL_MAX_SUPPORTED_HUFF_CODESIZE + 1]; DEFL_CLEAR_OBJ(num_codes);
		if (static_table)
		{
			for (i = 0; i < table_len; i++) num_codes[d->m_huff_code_sizes[table_num][i]]++;
		}
		else
		{
			defl_sym_freq syms0[DEFL_MAX_HUFF_SYMBOLS], syms1[DEFL_MAX_HUFF_SYMBOLS], * pSyms;
			int num_used_syms = 0;
			const uint16_t* pSym_count = &d->m_huff_count[table_num][0];
			for (i = 0; i < table_len; i++) if (pSym_count[i]) { syms0[num_used_syms].m_key = (uint16_t)pSym_count[i]; syms0[num_used_syms++].m_sym_index = (uint16_t)i; }

			pSyms = defl_radix_sort_syms(num_used_syms, syms0, syms1); defl_calculate_minimum_redundancy(pSyms, num_used_syms);

			for (i = 0; i < num_used_syms; i++) num_codes[pSyms[i].m_key]++;

			defl_huffman_enforce_max_code_size(num_codes, num_used_syms, code_size_limit);

			DEFL_CLEAR_OBJ(d->m_huff_code_sizes[table_num]); DEFL_CLEAR_OBJ(d->m_huff_codes[table_num]);
			for (i = 1, j = num_used_syms; i <= code_size_limit; i++)
				for (l = num_codes[i]; l > 0; l--) d->m_huff_code_sizes[table_num][pSyms[--j].m_sym_index] = (uint8_t)(i);
		}

		next_code[1] = 0; for (j = 0, i = 2; i <= code_size_limit; i++) next_code[i] = j = ((j + num_codes[i - 1]) << 1);

		for (i = 0; i < table_len; i++)
		{
			uint32_t rev_code = 0, code, code_size; if ((code_size = d->m_huff_code_sizes[table_num][i]) == 0) continue;
			code = next_code[code_size]++; for (l = code_size; l > 0; l--, code >>= 1) rev_code = (rev_code << 1) | (code & 1);
			d->m_huff_codes[table_num][i] = (uint16_t)rev_code;
		}
	}

#define DEFL_RLE_PREV_CODE_SIZE() { if (rle_repeat_count) { \
  if (rle_repeat_count < 3) { \
    d->m_huff_count[2][prev_code_size] = (uint16_t)(d->m_huff_count[2][prev_code_size] + rle_repeat_count); \
    while (rle_repeat_count--) packed_code_sizes[num_packed_code_sizes++] = prev_code_size; \
  } else { \
    d->m_huff_count[2][16] = (uint16_t)(d->m_huff_count[2][16] + 1); packed_code_sizes[num_packed_code_sizes++] = 16; packed_code_sizes[num_packed_code_sizes++] = (uint8_t)(rle_repeat_count - 3); \
} rle_repeat_count = 0; } }

#define DEFL_RLE_ZERO_CODE_SIZE() { if (rle_z_count) { \
  if (rle_z_count < 3) { \
    d->m_huff_count[2][0] = (uint16_t)(d->m_huff_count[2][0] + rle_z_count); while (rle_z_count--) packed_code_sizes[num_packed_code_sizes++] = 0; \
  } else if (rle_z_count <= 10) { \
    d->m_huff_count[2][17] = (uint16_t)(d->m_huff_count[2][17] + 1); packed_code_sizes[num_packed_code_sizes++] = 17; packed_code_sizes[num_packed_code_sizes++] = (uint8_t)(rle_z_count - 3); \
  } else { \
    d->m_huff_count[2][18] = (uint16_t)(d->m_huff_count[2][18] + 1); packed_code_sizes[num_packed_code_sizes++] = 18; packed_code_sizes[num_packed_code_sizes++] = (uint8_t)(rle_z_count - 11); \
} rle_z_count = 0; } }

	static uint8_t g_defl_packed_code_size_syms_swizzle[] = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };

#define DEFL_DYN_PUT_BITS(bb, ll) \
do { \
	uint32_t b = (bb), l = (ll); \
	assert((l) >= 1 && (l) <= 16); assert((b) < (1ULL << (l))); \
	bit_buf |= (((uint64_t)(b)) << bit_buf_size); bit_buf_size += (l); assert(bit_buf_size <= 64); \
	while (bit_buf_size >= 8) \
	{ \
		if ((dst_ofs + 1) > dst_buf_size) \
			return false; \
		*(uint8_t*)(pDst + dst_ofs) = (uint8_t)bit_buf; \
		dst_ofs++; \
		bit_buf >>= 8; \
		bit_buf_size -= 8; \
	} \
} while(0)

	static bool defl_start_dynamic_block(defl_huff* d, uint8_t* pDst, uint32_t& dst_ofs, uint32_t dst_buf_size, uint64_t& bit_buf, int& bit_buf_size)
	{
		int num_lit_codes, num_dist_codes, num_bit_lengths; uint32_t i, total_code_sizes_to_pack, num_packed_code_sizes, rle_z_count, rle_repeat_count, packed_code_sizes_index;
		uint8_t code_sizes_to_pack[DEFL_MAX_HUFF_SYMBOLS_0 + DEFL_MAX_HUFF_SYMBOLS_1], packed_code_sizes[DEFL_MAX_HUFF_SYMBOLS_0 + DEFL_MAX_HUFF_SYMBOLS_1], prev_code_size = 0xFF;

#if FPNG_TRAIN_HUFFMAN_TABLES
		assert(HUFF_COUNTS_SIZE == DEFL_MAX_HUFF_SYMBOLS_0);
		for (uint32_t i = 0; i < DEFL_MAX_HUFF_SYMBOLS_0; i++)
			g_huff_counts[i] += d->m_huff_count[0][i];
#endif

		d->m_huff_count[0][256] = 1;

		defl_optimize_huffman_table(d, 0, DEFL_MAX_HUFF_SYMBOLS_0, 12, FPNG_FALSE);
		defl_optimize_huffman_table(d, 1, DEFL_MAX_HUFF_SYMBOLS_1, 12, FPNG_FALSE);

		for (num_lit_codes = 286; num_lit_codes > 257; num_lit_codes--) if (d->m_huff_code_sizes[0][num_lit_codes - 1]) break;
		for (num_dist_codes = 30; num_dist_codes > 1; num_dist_codes--) if (d->m_huff_code_sizes[1][num_dist_codes - 1]) break;

		memcpy(code_sizes_to_pack, &d->m_huff_code_sizes[0][0], num_lit_codes);
		memcpy(code_sizes_to_pack + num_lit_codes, &d->m_huff_code_sizes[1][0], num_dist_codes);
		total_code_sizes_to_pack = num_lit_codes + num_dist_codes; num_packed_code_sizes = 0; rle_z_count = 0; rle_repeat_count = 0;

		memset(&d->m_huff_count[2][0], 0, sizeof(d->m_huff_count[2][0]) * DEFL_MAX_HUFF_SYMBOLS_2);
		for (i = 0; i < total_code_sizes_to_pack; i++)
		{
			uint8_t code_size = code_sizes_to_pack[i];
			if (!code_size)
			{
				DEFL_RLE_PREV_CODE_SIZE();
				if (++rle_z_count == 138) { DEFL_RLE_ZERO_CODE_SIZE(); }
			}
			else
			{
				DEFL_RLE_ZERO_CODE_SIZE();
				if (code_size != prev_code_size)
				{
					DEFL_RLE_PREV_CODE_SIZE();
					d->m_huff_count[2][code_size] = (uint16_t)(d->m_huff_count[2][code_size] + 1); packed_code_sizes[num_packed_code_sizes++] = code_size;
				}
				else if (++rle_repeat_count == 6)
				{
					DEFL_RLE_PREV_CODE_SIZE();
				}
			}
			prev_code_size = code_size;
		}
		if (rle_repeat_count) { DEFL_RLE_PREV_CODE_SIZE(); }
		else { DEFL_RLE_ZERO_CODE_SIZE(); }

		defl_optimize_huffman_table(d, 2, DEFL_MAX_HUFF_SYMBOLS_2, 7, FPNG_FALSE);

		// max of 2+5+5+4+18*3+(288+32)*7=2310 bits
		DEFL_DYN_PUT_BITS(2, 2);

		DEFL_DYN_PUT_BITS(num_lit_codes - 257, 5);
		DEFL_DYN_PUT_BITS(num_dist_codes - 1, 5);

		for (num_bit_lengths = 18; num_bit_lengths >= 0; num_bit_lengths--) if (d->m_huff_code_sizes[2][g_defl_packed_code_size_syms_swizzle[num_bit_lengths]]) break;
		num_bit_lengths = maximum<int>(4, (num_bit_lengths + 1)); DEFL_DYN_PUT_BITS(num_bit_lengths - 4, 4);
		for (i = 0; (int)i < num_bit_lengths; i++) DEFL_DYN_PUT_BITS(d->m_huff_code_sizes[2][g_defl_packed_code_size_syms_swizzle[i]], 3);

		for (packed_code_sizes_index = 0; packed_code_sizes_index < num_packed_code_sizes; )
		{
			uint32_t code = packed_code_sizes[packed_code_sizes_index++]; assert(code < DEFL_MAX_HUFF_SYMBOLS_2);
			DEFL_DYN_PUT_BITS(d->m_huff_codes[2][code], d->m_huff_code_sizes[2][code]);
			if (code >= 16) DEFL_DYN_PUT_BITS(packed_code_sizes[packed_code_sizes_index++], "\02\03\07"[code - 16]);
		}

		return true;
	}

	static uint32_t write_raw_block(const uint8_t* pSrc, uint32_t src_len, uint8_t* pDst, uint32_t dst_buf_size)
	{
		if (dst_buf_size < 2)
			return 0;

		pDst[0] = 0x78;
		pDst[1] = 0x01;

		uint32_t dst_ofs = 2;

		uint32_t src_ofs = 0;
		while (src_ofs < src_len)
		{
			const uint32_t src_remaining = src_len - src_ofs;
			const uint32_t block_size = minimum<uint32_t>(UINT16_MAX, src_remaining);
			const bool final_block = (block_size == src_remaining);

			if ((dst_ofs + 5 + block_size) > dst_buf_size)
				return 0;

			pDst[dst_ofs + 0] = final_block ? 1 : 0;

			pDst[dst_ofs + 1] = block_size & 0xFF;
			pDst[dst_ofs + 2] = (block_size >> 8) & 0xFF;

			pDst[dst_ofs + 3] = (~block_size) & 0xFF;
			pDst[dst_ofs + 4] = ((~block_size) >> 8) & 0xFF;

			memcpy(pDst + dst_ofs + 5, pSrc + src_ofs, block_size);

			src_ofs += block_size;
			dst_ofs += 5 + block_size;
		}

		uint32_t src_adler32 = fpng_adler32(pSrc, src_len, FPNG_ADLER32_INIT);

		for (uint32_t i = 0; i < 4; i++)
		{
			if (dst_ofs + 1 > dst_buf_size)
				return 0;

			pDst[dst_ofs] = (uint8_t)(src_adler32 >> 24);
			dst_ofs++;

			src_adler32 <<= 8;
		}

		return dst_ofs;
	}

	static void adjust_freq32(uint32_t num_freq, uint32_t* pFreq, uint16_t* pFreq16)
	{
		uint32_t total_freq = 0;
		for (uint32_t i = 0; i < num_freq; i++)
			total_freq += pFreq[i];

		if (!total_freq)
		{
			memset(pFreq16, 0, num_freq * sizeof(uint16_t));
			return;
		}

		uint32_t total_freq16 = 0;
		for (uint32_t i = 0; i < num_freq; i++)
		{
			uint64_t f = pFreq[i];
			if (!f)
			{
				pFreq16[i] = 0;
				continue;
			}

			pFreq16[i] = (uint16_t)maximum<uint32_t>(1, (uint32_t)((f * UINT16_MAX) / total_freq));

			total_freq16 += pFreq16[i];
		}

		while (total_freq16 > UINT16_MAX)
		{
			total_freq16 = 0;
			for (uint32_t i = 0; i < num_freq; i++)
			{
				if (pFreq[i])
				{
					pFreq[i] = maximum<uint32_t>(1, pFreq[i] >> 1);
					total_freq16 += pFreq[i];
				}
			}
		}
	}

#if FPNG_TRAIN_HUFFMAN_TABLES
	bool create_dynamic_block_prefix(uint64_t* pFreq, uint32_t num_chans, std::vector<uint8_t>& prefix, uint64_t& bit_buf, int &bit_buf_size, uint32_t* pCodes, uint8_t* pCodesizes)
	{
		assert((num_chans == 3) || (num_chans == 4));
		assert(HUFF_COUNTS_SIZE == DEFL_MAX_HUFF_SYMBOLS_0); // must be equal
				
		defl_huff dh;
		memset(&dh, 0, sizeof(dh));

		uint32_t lit_freq[DEFL_MAX_HUFF_SYMBOLS_0];
		
		uint32_t shift_len = 0;
		for (; ; )
		{
			uint32_t i;
			for (i = 0; i < DEFL_MAX_HUFF_SYMBOLS_0; i++)
			{
				uint64_t f = pFreq[i];
				if (f)
					f = maximum<uint64_t>(1U, f >> shift_len);

				if (f > UINT32_MAX)
					break;

				lit_freq[i] = (uint32_t)pFreq[i];
			}

			if (i == DEFL_MAX_HUFF_SYMBOLS_0)
				break;
			
			shift_len++;
		}
				
		// Ensure all valid Deflate literal/EOB/length syms are non-zero, so anything can be coded.
		for (uint32_t i = 0; i <= 256; i++)
		{
			if (!lit_freq[i])
				lit_freq[i] = 1;
		}

		for (uint32_t len = num_chans; len <= DEFL_MAX_MATCH_LEN; len += num_chans)
		{
			uint32_t sym = g_defl_len_sym[len - 3];
			if (!lit_freq[sym])
				lit_freq[sym] = 1;
		}

		adjust_freq32(DEFL_MAX_HUFF_SYMBOLS_0, lit_freq, &dh.m_huff_count[0][0]);
		
		const uint32_t dist_sym = g_defl_small_dist_sym[num_chans - 1];
		dh.m_huff_count[1][dist_sym] = 1;
		dh.m_huff_count[1][dist_sym + 1] = 1; // to workaround a bug in wuffs decoder
			
		prefix.resize(4096);
		uint8_t* pDst = prefix.data();
		uint32_t dst_buf_size = (uint32_t)prefix.size();

		uint32_t dst_ofs = 0;

		// zlib header
		PUT_BITS(0x78, 8);
		PUT_BITS(0x01, 8);

		// write BFINAL bit
		PUT_BITS(1, 1);
				
		if (!defl_start_dynamic_block(&dh, pDst, dst_ofs, dst_buf_size, bit_buf, bit_buf_size))
			return false;

		prefix.resize(dst_ofs);

		for (uint32_t i = 0; i < DEFL_MAX_HUFF_SYMBOLS_0; i++)
		{
			pCodes[i] = dh.m_huff_codes[0][i];
			pCodesizes[i] = dh.m_huff_code_sizes[0][i];
		}

		return true;
	}
#endif

	static uint32_t pixel_deflate_dyn_3_rle(
		const uint8_t* pImg, uint32_t w, uint32_t h,
		uint8_t* pDst, uint32_t dst_buf_size)
	{
		const uint32_t bpl = 1 + w * 3;

		uint64_t bit_buf = 0;
		int bit_buf_size = 0;

		uint32_t dst_ofs = 0;

		// zlib header
		PUT_BITS(0x78, 8);
		PUT_BITS(0x01, 8);

		// write BFINAL bit
		PUT_BITS(1, 1);

		std::vector<uint32_t> codes((w + 1) * h);
		uint32_t* pDst_codes = codes.data();

		uint32_t lit_freq[DEFL_MAX_HUFF_SYMBOLS_0];
		memset(lit_freq, 0, sizeof(lit_freq));
		
		const uint8_t* pSrc = pImg;
		uint32_t src_ofs = 0;

		uint32_t src_adler32 = fpng_adler32(pImg, bpl * h, FPNG_ADLER32_INIT);

		const uint32_t dist_sym = g_defl_small_dist_sym[3 - 1];
				
		for (uint32_t y = 0; y < h; y++)
		{
			const uint32_t end_src_ofs = src_ofs + bpl;

			const uint32_t filter_lit = pSrc[src_ofs++];
			*pDst_codes++ = 1 | (filter_lit << 8);
			lit_freq[filter_lit]++;

			uint32_t prev_lits;

			{
				uint32_t lits = READ_RGB_PIXEL(pSrc + src_ofs);

				*pDst_codes++ = lits << 8;

				lit_freq[lits & 0xFF]++;
				lit_freq[(lits >> 8) & 0xFF]++;
				lit_freq[lits >> 16]++;

				src_ofs += 3;

				prev_lits = lits;
			}

			while (src_ofs < end_src_ofs)
			{
				uint32_t lits = READ_RGB_PIXEL(pSrc + src_ofs);

				if (lits == prev_lits)
				{
					uint32_t match_len = 3;
					uint32_t max_match_len = minimum<int>(255, (int)(end_src_ofs - src_ofs));

					while (match_len < max_match_len)
					{
						if (READ_RGB_PIXEL(pSrc + src_ofs + match_len) != lits)
							break;
						match_len += 3;
					}
										
					*pDst_codes++ = match_len - 1;

					uint32_t adj_match_len = match_len - 3;

					lit_freq[g_defl_len_sym[adj_match_len]]++;
					
					src_ofs += match_len;
				}
				else
				{
					*pDst_codes++ = lits << 8;

					lit_freq[lits & 0xFF]++;
					lit_freq[(lits >> 8) & 0xFF]++;
					lit_freq[lits >> 16]++;

					prev_lits = lits;

					src_ofs += 3;
				}

			} // while (src_ofs < end_src_ofs)

		} // y

		assert(src_ofs == h * bpl);
		const uint32_t total_codes = (uint32_t)(pDst_codes - codes.data());
		assert(total_codes <= codes.size());
								
		defl_huff dh;
		
		lit_freq[256] = 1;

		adjust_freq32(DEFL_MAX_HUFF_SYMBOLS_0, lit_freq, &dh.m_huff_count[0][0]);

		memset(&dh.m_huff_count[1][0], 0, sizeof(dh.m_huff_count[1][0]) * DEFL_MAX_HUFF_SYMBOLS_1);
		dh.m_huff_count[1][dist_sym] = 1;
		dh.m_huff_count[1][dist_sym + 1] = 1; // to workaround a bug in wuffs decoder

		if (!defl_start_dynamic_block(&dh, pDst, dst_ofs, dst_buf_size, bit_buf, bit_buf_size))
			return 0;

		assert(bit_buf_size <= 7);
		assert(dh.m_huff_codes[1][dist_sym] == 0 && dh.m_huff_code_sizes[1][dist_sym] == 1);
				
		for (uint32_t i = 0; i < total_codes; i++)
		{
			uint32_t c = codes[i];

			uint32_t c_type = c & 0xFF;
			if (c_type == 0)
			{
				uint32_t lits = c >> 8;

				PUT_BITS_CZ(dh.m_huff_codes[0][lits & 0xFF], dh.m_huff_code_sizes[0][lits & 0xFF]);
				lits >>= 8;

				PUT_BITS_CZ(dh.m_huff_codes[0][lits & 0xFF], dh.m_huff_code_sizes[0][lits & 0xFF]);
				lits >>= 8;

				PUT_BITS_CZ(dh.m_huff_codes[0][lits], dh.m_huff_code_sizes[0][lits]);
			}
			else if (c_type == 1)
			{
				uint32_t lit = c >> 8;
				PUT_BITS_CZ(dh.m_huff_codes[0][lit], dh.m_huff_code_sizes[0][lit]);
			}
			else
			{
				uint32_t match_len = c_type + 1;

				uint32_t adj_match_len = match_len - 3;
				
				PUT_BITS_CZ(dh.m_huff_codes[0][g_defl_len_sym[adj_match_len]], dh.m_huff_code_sizes[0][g_defl_len_sym[adj_match_len]]);
				PUT_BITS(adj_match_len & g_bitmasks[g_defl_len_extra[adj_match_len]], g_defl_len_extra[adj_match_len] + 1); // up to 6 bits, +1 for the match distance Huff code which is always 0

				// no need to write the distance code, it's always 0
				//PUT_BITS_CZ(dh.m_huff_codes[1][dist_sym], dh.m_huff_code_sizes[1][dist_sym]);
			}

			// up to 55 bits
			PUT_BITS_FLUSH;
		}

		PUT_BITS_CZ(dh.m_huff_codes[0][256], dh.m_huff_code_sizes[0][256]);

		PUT_BITS_FORCE_FLUSH;

		// Write zlib adler32
		for (uint32_t i = 0; i < 4; i++)
		{
			if ((dst_ofs + 1) > dst_buf_size)
				return 0;
			*(uint8_t*)(pDst + dst_ofs) = (uint8_t)(src_adler32 >> 24);
			dst_ofs++;

			src_adler32 <<= 8;
		}

		return dst_ofs;
	}

	static uint32_t pixel_deflate_dyn_3_rle_one_pass(
		const uint8_t* pImg, uint32_t w, uint32_t h,
		uint8_t* pDst, uint32_t dst_buf_size)
	{
		const uint32_t bpl = 1 + w * 3;

		if (dst_buf_size < sizeof(g_dyn_huff_3))
			return false;
		memcpy(pDst, g_dyn_huff_3, sizeof(g_dyn_huff_3));
		uint32_t dst_ofs = sizeof(g_dyn_huff_3);

		uint64_t bit_buf = DYN_HUFF_3_BITBUF;
		int bit_buf_size = DYN_HUFF_3_BITBUF_SIZE;

		const uint8_t* pSrc = pImg;
		uint32_t src_ofs = 0;

		uint32_t src_adler32 = fpng_adler32(pImg, bpl * h, FPNG_ADLER32_INIT);

		for (uint32_t y = 0; y < h; y++)
		{
			const uint32_t end_src_ofs = src_ofs + bpl;

			const uint32_t filter_lit = pSrc[src_ofs++];
			PUT_BITS_CZ(g_dyn_huff_3_codes[filter_lit].m_code, g_dyn_huff_3_codes[filter_lit].m_code_size);

			uint32_t prev_lits;

			{
				uint32_t lits = READ_RGB_PIXEL(pSrc + src_ofs);

				PUT_BITS_CZ(g_dyn_huff_3_codes[lits & 0xFF].m_code, g_dyn_huff_3_codes[lits & 0xFF].m_code_size);
				PUT_BITS_CZ(g_dyn_huff_3_codes[(lits >> 8) & 0xFF].m_code, g_dyn_huff_3_codes[(lits >> 8) & 0xFF].m_code_size);
				PUT_BITS_CZ(g_dyn_huff_3_codes[(lits >> 16)].m_code, g_dyn_huff_3_codes[(lits >> 16)].m_code_size);

				src_ofs += 3;
			
				prev_lits = lits;
			}

			PUT_BITS_FLUSH;

			while (src_ofs < end_src_ofs)
			{
				uint32_t lits = READ_RGB_PIXEL(pSrc + src_ofs);

				if (lits == prev_lits)
				{
					uint32_t match_len = 3;
					uint32_t max_match_len = minimum<int>(255, (int)(end_src_ofs - src_ofs));

					while (match_len < max_match_len)
					{
						if (READ_RGB_PIXEL(pSrc + src_ofs + match_len) != lits)
							break;
						match_len += 3;
					}
										
					uint32_t adj_match_len = match_len - 3;

					PUT_BITS_CZ(g_dyn_huff_3_codes[g_defl_len_sym[adj_match_len]].m_code, g_dyn_huff_3_codes[g_defl_len_sym[adj_match_len]].m_code_size);
					PUT_BITS(adj_match_len & g_bitmasks[g_defl_len_extra[adj_match_len]], g_defl_len_extra[adj_match_len] + 1); // up to 6 bits, +1 for the match distance Huff code which is always 0

					src_ofs += match_len;
				}
				else
				{
					PUT_BITS_CZ(g_dyn_huff_3_codes[lits & 0xFF].m_code, g_dyn_huff_3_codes[lits & 0xFF].m_code_size);
					PUT_BITS_CZ(g_dyn_huff_3_codes[(lits >> 8) & 0xFF].m_code, g_dyn_huff_3_codes[(lits >> 8) & 0xFF].m_code_size);
					PUT_BITS_CZ(g_dyn_huff_3_codes[(lits >> 16)].m_code, g_dyn_huff_3_codes[(lits >> 16)].m_code_size);
					
					prev_lits = lits;

					src_ofs += 3;
				}

				PUT_BITS_FLUSH;

			} // while (src_ofs < end_src_ofs)

		} // y

		assert(src_ofs == h * bpl);
		
		assert(bit_buf_size <= 7);

		PUT_BITS_CZ(g_dyn_huff_3_codes[256].m_code, g_dyn_huff_3_codes[256].m_code_size);

		PUT_BITS_FORCE_FLUSH;

		// Write zlib adler32
		for (uint32_t i = 0; i < 4; i++)
		{
			if ((dst_ofs + 1) > dst_buf_size)
				return 0;
			*(uint8_t*)(pDst + dst_ofs) = (uint8_t)(src_adler32 >> 24);
			dst_ofs++;

			src_adler32 <<= 8;
		}

		return dst_ofs;
	}

	static uint32_t pixel_deflate_dyn_4_rle(
		const uint8_t* pImg, uint32_t w, uint32_t h,
		uint8_t* pDst, uint32_t dst_buf_size)
	{
		const uint32_t bpl = 1 + w * 4;

		uint64_t bit_buf = 0;
		int bit_buf_size = 0;

		uint32_t dst_ofs = 0;

		// zlib header
		PUT_BITS(0x78, 8);
		PUT_BITS(0x01, 8);

		// write BFINAL bit
		PUT_BITS(1, 1);

		std::vector<uint64_t> codes;
		codes.resize((w + 1) * h);
		uint64_t* pDst_codes = codes.data();

		uint32_t lit_freq[DEFL_MAX_HUFF_SYMBOLS_0];
		memset(lit_freq, 0, sizeof(lit_freq));

		const uint8_t* pSrc = pImg;
		uint32_t src_ofs = 0;

		uint32_t src_adler32 = fpng_adler32(pImg, bpl * h, FPNG_ADLER32_INIT);

		const uint32_t dist_sym = g_defl_small_dist_sym[4 - 1];

		for (uint32_t y = 0; y < h; y++)
		{
			const uint32_t end_src_ofs = src_ofs + bpl;

			const uint32_t filter_lit = pSrc[src_ofs++];
			*pDst_codes++ = 1 | (filter_lit << 8);
			lit_freq[filter_lit]++;

			uint32_t prev_lits;
			{
				uint32_t lits = READ_LE32(pSrc + src_ofs);

				*pDst_codes++ = (uint64_t)lits << 8;

				lit_freq[lits & 0xFF]++;
				lit_freq[(lits >> 8) & 0xFF]++;
				lit_freq[(lits >> 16) & 0xFF]++;
				lit_freq[lits >> 24]++;

				src_ofs += 4;
				
				prev_lits = lits;
			}

			while (src_ofs < end_src_ofs)
			{
				uint32_t lits = READ_LE32(pSrc + src_ofs);

				if (lits == prev_lits)
				{
					uint32_t match_len = 4;
					uint32_t max_match_len = minimum<int>(252, (int)(end_src_ofs - src_ofs));

					while (match_len < max_match_len)
					{
						if (READ_LE32(pSrc + src_ofs + match_len) != lits)
							break;
						match_len += 4;
					}
										
					*pDst_codes++ = match_len - 1;

					uint32_t adj_match_len = match_len - 3;

					lit_freq[g_defl_len_sym[adj_match_len]]++;
					
					src_ofs += match_len;
				}
				else
				{
					*pDst_codes++ = (uint64_t)lits << 8;

					lit_freq[lits & 0xFF]++;
					lit_freq[(lits >> 8) & 0xFF]++;
					lit_freq[(lits >> 16) & 0xFF]++;
					lit_freq[lits >> 24]++;
					
					prev_lits = lits;

					src_ofs += 4;
				}

			} // while (src_ofs < end_src_ofs)

		} // y

		assert(src_ofs == h * bpl);
		const uint32_t total_codes = (uint32_t)(pDst_codes - codes.data());
		assert(total_codes <= codes.size());
						
		defl_huff dh;
		
		lit_freq[256] = 1;

		adjust_freq32(DEFL_MAX_HUFF_SYMBOLS_0, lit_freq, &dh.m_huff_count[0][0]);
		
		memset(&dh.m_huff_count[1][0], 0, sizeof(dh.m_huff_count[1][0]) * DEFL_MAX_HUFF_SYMBOLS_1);
		dh.m_huff_count[1][dist_sym] = 1;
		dh.m_huff_count[1][dist_sym + 1] = 1; // to workaround a bug in wuffs decoder

		if (!defl_start_dynamic_block(&dh, pDst, dst_ofs, dst_buf_size, bit_buf, bit_buf_size))
			return 0;

		assert(bit_buf_size <= 7);
		assert(dh.m_huff_codes[1][dist_sym] == 0 && dh.m_huff_code_sizes[1][dist_sym] == 1);

		for (uint32_t i = 0; i < total_codes; i++)
		{
			uint64_t c = codes[i];

			uint32_t c_type = (uint32_t)(c & 0xFF);
			if (c_type == 0)
			{
				uint32_t lits = (uint32_t)(c >> 8);

				PUT_BITS_CZ(dh.m_huff_codes[0][lits & 0xFF], dh.m_huff_code_sizes[0][lits & 0xFF]);
				lits >>= 8;

				PUT_BITS_CZ(dh.m_huff_codes[0][lits & 0xFF], dh.m_huff_code_sizes[0][lits & 0xFF]);
				lits >>= 8;

				PUT_BITS_CZ(dh.m_huff_codes[0][lits & 0xFF], dh.m_huff_code_sizes[0][lits & 0xFF]);
				lits >>= 8;

				if (bit_buf_size >= 49)
				{
					PUT_BITS_FLUSH;
				}

				PUT_BITS_CZ(dh.m_huff_codes[0][lits], dh.m_huff_code_sizes[0][lits]);
			}
			else if (c_type == 1)
			{
				uint32_t lit = (uint32_t)(c >> 8);
				PUT_BITS_CZ(dh.m_huff_codes[0][lit], dh.m_huff_code_sizes[0][lit]);
			}
			else
			{
				uint32_t match_len = c_type + 1;

				uint32_t adj_match_len = match_len - 3;
				
				PUT_BITS_CZ(dh.m_huff_codes[0][g_defl_len_sym[adj_match_len]], dh.m_huff_code_sizes[0][g_defl_len_sym[adj_match_len]]);
				PUT_BITS(adj_match_len & g_bitmasks[g_defl_len_extra[adj_match_len]], g_defl_len_extra[adj_match_len] + 1); // up to 6 bits, +1 for the match distance Huff code which is always 0

				// no need to write the distance code, it's always 0
			}

			// up to 55 bits
			PUT_BITS_FLUSH;
		}

		PUT_BITS_CZ(dh.m_huff_codes[0][256], dh.m_huff_code_sizes[0][256]);

		PUT_BITS_FORCE_FLUSH;

		// Write zlib adler32
		for (uint32_t i = 0; i < 4; i++)
		{
			if ((dst_ofs + 1) > dst_buf_size)
				return 0;
			*(uint8_t*)(pDst + dst_ofs) = (uint8_t)(src_adler32 >> 24);
			dst_ofs++;

			src_adler32 <<= 8;
		}

		return dst_ofs;
	}

	static uint32_t pixel_deflate_dyn_4_rle_one_pass(
		const uint8_t* pImg, uint32_t w, uint32_t h,
		uint8_t* pDst, uint32_t dst_buf_size)
	{
		const uint32_t bpl = 1 + w * 4;

		if (dst_buf_size < sizeof(g_dyn_huff_4))
			return false;
		memcpy(pDst, g_dyn_huff_4, sizeof(g_dyn_huff_4));
		uint32_t dst_ofs = sizeof(g_dyn_huff_4);

		uint64_t bit_buf = DYN_HUFF_4_BITBUF;
		int bit_buf_size = DYN_HUFF_4_BITBUF_SIZE;

		const uint8_t* pSrc = pImg;
		uint32_t src_ofs = 0;

		uint32_t src_adler32 = fpng_adler32(pImg, bpl * h, FPNG_ADLER32_INIT);

		for (uint32_t y = 0; y < h; y++)
		{
			const uint32_t end_src_ofs = src_ofs + bpl;

			const uint32_t filter_lit = pSrc[src_ofs++];
			PUT_BITS_CZ(g_dyn_huff_4_codes[filter_lit].m_code, g_dyn_huff_4_codes[filter_lit].m_code_size);

			PUT_BITS_FLUSH;

			uint32_t prev_lits;
			{
				uint32_t lits = READ_LE32(pSrc + src_ofs);

				PUT_BITS_CZ(g_dyn_huff_4_codes[lits & 0xFF].m_code, g_dyn_huff_4_codes[lits & 0xFF].m_code_size);
				PUT_BITS_CZ(g_dyn_huff_4_codes[(lits >> 8) & 0xFF].m_code, g_dyn_huff_4_codes[(lits >> 8) & 0xFF].m_code_size);
				PUT_BITS_CZ(g_dyn_huff_4_codes[(lits >> 16) & 0xFF].m_code, g_dyn_huff_4_codes[(lits >> 16) & 0xFF].m_code_size);

				if (bit_buf_size >= 49)
				{
					PUT_BITS_FLUSH;
				}
				
				PUT_BITS_CZ(g_dyn_huff_4_codes[(lits >> 24)].m_code, g_dyn_huff_4_codes[(lits >> 24)].m_code_size);

				src_ofs += 4;
				
				prev_lits = lits;
			}

			PUT_BITS_FLUSH;

			while (src_ofs < end_src_ofs)
			{
				uint32_t lits = READ_LE32(pSrc + src_ofs);
								
				if (lits == prev_lits)
				{
					uint32_t match_len = 4;
					uint32_t max_match_len = minimum<int>(252, (int)(end_src_ofs - src_ofs));

					while (match_len < max_match_len)
					{
						if (READ_LE32(pSrc + src_ofs + match_len) != lits)
							break;
						match_len += 4;
					}

					uint32_t adj_match_len = match_len - 3;

					const uint32_t match_code_bits = g_dyn_huff_4_codes[g_defl_len_sym[adj_match_len]].m_code_size;
					const uint32_t len_extra_bits = g_defl_len_extra[adj_match_len];

					if (match_len == 4)
					{
						// This check is optional - see if just encoding 4 literals would be cheaper than using a short match.
						uint32_t lit_bits = g_dyn_huff_4_codes[lits & 0xFF].m_code_size + g_dyn_huff_4_codes[(lits >> 8) & 0xFF].m_code_size + 
							g_dyn_huff_4_codes[(lits >> 16) & 0xFF].m_code_size + g_dyn_huff_4_codes[(lits >> 24)].m_code_size;
						
						if ((match_code_bits + len_extra_bits + 1) > lit_bits)
							goto do_literals;
					}

					PUT_BITS_CZ(g_dyn_huff_4_codes[g_defl_len_sym[adj_match_len]].m_code, match_code_bits);
					PUT_BITS(adj_match_len & g_bitmasks[g_defl_len_extra[adj_match_len]], len_extra_bits + 1); // up to 6 bits, +1 for the match distance Huff code which is always 0

					src_ofs += match_len;
				}
				else
				{
do_literals:
					PUT_BITS_CZ(g_dyn_huff_4_codes[lits & 0xFF].m_code, g_dyn_huff_4_codes[lits & 0xFF].m_code_size);
					PUT_BITS_CZ(g_dyn_huff_4_codes[(lits >> 8) & 0xFF].m_code, g_dyn_huff_4_codes[(lits >> 8) & 0xFF].m_code_size);
					PUT_BITS_CZ(g_dyn_huff_4_codes[(lits >> 16) & 0xFF].m_code, g_dyn_huff_4_codes[(lits >> 16) & 0xFF].m_code_size);

					if (bit_buf_size >= 49)
					{
						PUT_BITS_FLUSH;
					}

					PUT_BITS_CZ(g_dyn_huff_4_codes[(lits >> 24)].m_code, g_dyn_huff_4_codes[(lits >> 24)].m_code_size);

					src_ofs += 4;
					
					prev_lits = lits;
				}

				PUT_BITS_FLUSH;

			} // while (src_ofs < end_src_ofs)

		} // y

		assert(src_ofs == h * bpl);

		assert(bit_buf_size <= 7);

		PUT_BITS_CZ(g_dyn_huff_4_codes[256].m_code, g_dyn_huff_4_codes[256].m_code_size);

		PUT_BITS_FORCE_FLUSH;

		// Write zlib adler32
		for (uint32_t i = 0; i < 4; i++)
		{
			if ((dst_ofs + 1) > dst_buf_size)
				return 0;
			*(uint8_t*)(pDst + dst_ofs) = (uint8_t)(src_adler32 >> 24);
			dst_ofs++;

			src_adler32 <<= 8;
		}

		return dst_ofs;
	}

	static void vector_append(std::vector<uint8_t>& buf, const void* pData, size_t len)
	{
		if (len)
		{
			size_t l = buf.size();
			buf.resize(l + len);
			memcpy(buf.data() + l, pData, len);
		}
	}
		
	static void apply_filter(uint32_t filter, int w, int h, uint32_t num_chans, uint32_t bpl, const uint8_t* pSrc, const uint8_t* pPrev_src, uint8_t* pDst)
	{
		(void)h;

		switch (filter)
		{
		case 0:
		{
			*pDst++ = 0;

			memcpy(pDst, pSrc, bpl);
			break;
		}
		case 2:
		{
			assert(pPrev_src);

			// Previous scanline
			*pDst++ = 2;

#if FPNG_X86_OR_X64_CPU && !FPNG_NO_SSE
			if (g_cpu_info.can_use_sse41())
			{
				uint32_t bytes_to_process = w * num_chans, ofs = 0;
				for (; bytes_to_process >= 16; bytes_to_process -= 16, ofs += 16)
					_mm_storeu_si128((__m128i*)(pDst + ofs), _mm_sub_epi8(_mm_loadu_si128((const __m128i*)(pSrc + ofs)), _mm_loadu_si128((const __m128i*)(pPrev_src + ofs))));

				for (; bytes_to_process; bytes_to_process--, ofs++)
					pDst[ofs] = (uint8_t)(pSrc[ofs] - pPrev_src[ofs]);
			}
			else
#endif
			{
				if (num_chans == 3)
				{
					for (uint32_t x = 0; x < (uint32_t)w; x++)
					{
						pDst[0] = (uint8_t)(pSrc[0] - pPrev_src[0]);
						pDst[1] = (uint8_t)(pSrc[1] - pPrev_src[1]);
						pDst[2] = (uint8_t)(pSrc[2] - pPrev_src[2]);

						pSrc += 3;
						pPrev_src += 3;
						pDst += 3;
					}
				}
				else
				{
					for (uint32_t x = 0; x < (uint32_t)w; x++)
					{
						pDst[0] = (uint8_t)(pSrc[0] - pPrev_src[0]);
						pDst[1] = (uint8_t)(pSrc[1] - pPrev_src[1]);
						pDst[2] = (uint8_t)(pSrc[2] - pPrev_src[2]);
						pDst[3] = (uint8_t)(pSrc[3] - pPrev_src[3]);

						pSrc += 4;
						pPrev_src += 4;
						pDst += 4;
					}
				}
			}

			break;
		}
		default:
			assert(0);
			break;
		}
	}

	bool fpng_encode_image_to_memory(const void* pImage, uint32_t w, uint32_t h, uint32_t num_chans, std::vector<uint8_t>& out_buf, uint32_t flags)
	{
		if (!endian_check())
		{
			assert(0);
			return false;
		}

		if ((w < 1) || (h < 1) || (w * (uint64_t)h > UINT32_MAX) || (w > FPNG_MAX_SUPPORTED_DIM) || (h > FPNG_MAX_SUPPORTED_DIM))
		{
			assert(0);
			return false;
		}

		if ((num_chans != 3) && (num_chans != 4))
		{
			assert(0);
			return false;
		}

		int i, bpl = w * num_chans;
		uint32_t y;

		std::vector<uint8_t> temp_buf;
		temp_buf.resize((bpl + 1) * h + 7);
		uint32_t temp_buf_ofs = 0;

		for (y = 0; y < h; ++y)
		{
			const uint8_t* pSrc = (uint8_t*)pImage + y * bpl;
			const uint8_t* pPrev_src = y ? ((uint8_t*)pImage + (y - 1) * bpl) : nullptr;

			uint8_t* pDst = &temp_buf[temp_buf_ofs];

			apply_filter(y ? 2 : 0, w, h, num_chans, bpl, pSrc, pPrev_src, pDst);

			temp_buf_ofs += 1 + bpl;
		}

		const uint32_t PNG_HEADER_SIZE = 58;
				
		uint32_t out_ofs = PNG_HEADER_SIZE;
				
		out_buf.resize((out_ofs + (bpl + 1) * h + 7) & ~7);

		uint32_t defl_size = 0;
		if ((flags & FPNG_FORCE_UNCOMPRESSED) == 0)
		{
			if (num_chans == 3)
			{
				if (flags & FPNG_ENCODE_SLOWER)
					defl_size = pixel_deflate_dyn_3_rle(temp_buf.data(), w, h, &out_buf[out_ofs], (uint32_t)out_buf.size() - out_ofs);
				else
					defl_size = pixel_deflate_dyn_3_rle_one_pass(temp_buf.data(), w, h, &out_buf[out_ofs], (uint32_t)out_buf.size() - out_ofs);
			}
			else
			{
				if (flags & FPNG_ENCODE_SLOWER)
					defl_size = pixel_deflate_dyn_4_rle(temp_buf.data(), w, h, &out_buf[out_ofs], (uint32_t)out_buf.size() - out_ofs);
				else
					defl_size = pixel_deflate_dyn_4_rle_one_pass(temp_buf.data(), w, h, &out_buf[out_ofs], (uint32_t)out_buf.size() - out_ofs);
			}
		}

		uint32_t zlib_size = defl_size;
		
		if (!defl_size)
		{
			// Dynamic block failed to compress - fall back to uncompressed blocks, filter 0.

			temp_buf_ofs = 0;

			for (y = 0; y < h; ++y)
			{
				const uint8_t* pSrc = (uint8_t*)pImage + y * bpl;

				uint8_t* pDst = &temp_buf[temp_buf_ofs];

				apply_filter(0, w, h, num_chans, bpl, pSrc, nullptr, pDst);

				temp_buf_ofs += 1 + bpl;
			}

			assert(temp_buf_ofs <= temp_buf.size());
						
			out_buf.resize(out_ofs + 6 + temp_buf_ofs + ((temp_buf_ofs + 65534) / 65535) * 5);

			uint32_t raw_size = write_raw_block(temp_buf.data(), (uint32_t)temp_buf_ofs, out_buf.data() + out_ofs, (uint32_t)out_buf.size() - out_ofs);
			if (!raw_size)
			{
				// Somehow we miscomputed the size of the output buffer.
				assert(0);
				return false;
			}

			zlib_size = raw_size;
		}
		
		assert((out_ofs + zlib_size) <= out_buf.size());

		out_buf.resize(out_ofs + zlib_size);

		const uint32_t idat_len = (uint32_t)out_buf.size() - PNG_HEADER_SIZE;

		// Write real PNG header, fdEC chunk, and the beginning of the IDAT chunk
		{
			static const uint8_t s_color_type[] = { 0x00, 0x00, 0x04, 0x02, 0x06 };

			uint8_t pnghdr[58] = { 
				0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a,   // PNG sig
				0x00,0x00,0x00,0x0d, 'I','H','D','R',  // IHDR chunk len, type
			    0,0,(uint8_t)(w >> 8),(uint8_t)w, // width
				0,0,(uint8_t)(h >> 8),(uint8_t)h, // height
				8,   //bit_depth
				s_color_type[num_chans], // color_type
				0, // compression
				0, // filter
				0, // interlace
				0, 0, 0, 0, // IHDR crc32
				0, 0, 0, 5, 'f', 'd', 'E', 'C', 82, 36, 147, 227, FPNG_FDEC_VERSION,   0xE5, 0xAB, 0x62, 0x99, // our custom private, ancillary, do not copy, fdEC chunk
			  (uint8_t)(idat_len >> 24),(uint8_t)(idat_len >> 16),(uint8_t)(idat_len >> 8),(uint8_t)idat_len, 'I','D','A','T' // IDATA chunk len, type
			}; 

			// Compute IHDR CRC32
			uint32_t c = (uint32_t)fpng_crc32(pnghdr + 12, 17, FPNG_CRC32_INIT);
			for (i = 0; i < 4; ++i, c <<= 8)
				((uint8_t*)(pnghdr + 29))[i] = (uint8_t)(c >> 24);

			memcpy(out_buf.data(), pnghdr, PNG_HEADER_SIZE);
		}

		// Write IDAT chunk's CRC32 and a 0 length IEND chunk
		vector_append(out_buf, "\0\0\0\0\0\0\0\0\x49\x45\x4e\x44\xae\x42\x60\x82", 16); // IDAT CRC32, followed by the IEND chunk

		// Compute IDAT crc32
		uint32_t c = (uint32_t)fpng_crc32(out_buf.data() + PNG_HEADER_SIZE - 4, idat_len + 4, FPNG_CRC32_INIT);
		
		for (i = 0; i < 4; ++i, c <<= 8)
			(out_buf.data() + out_buf.size() - 16)[i] = (uint8_t)(c >> 24);
				
		return true;
	}

#ifndef FPNG_NO_STDIO
	bool fpng_encode_image_to_file(const char* pFilename, const void* pImage, uint32_t w, uint32_t h, uint32_t num_chans, uint32_t flags)
	{
		std::vector<uint8_t> out_buf;
		if (!fpng_encode_image_to_memory(pImage, w, h, num_chans, out_buf, flags))
			return false;

		FILE* pFile = nullptr;
#ifdef _MSC_VER
		fopen_s(&pFile, pFilename, "wb");
#else
		pFile = fopen(pFilename, "wb");
#endif
		if (!pFile)
			return false;

		if (fwrite(out_buf.data(), 1, out_buf.size(), pFile) != out_buf.size())
		{
			fclose(pFile);
			return false;
		}

		return (fclose(pFile) != EOF);
	}
#endif

	// Decompression

	const uint32_t FPNG_DECODER_TABLE_BITS = 12;
	const uint32_t FPNG_DECODER_TABLE_SIZE = 1 << FPNG_DECODER_TABLE_BITS;

	static bool build_decoder_table(uint32_t num_syms, uint8_t* pCode_sizes, uint32_t* pTable)
	{
		uint32_t num_codes[16];

		memset(num_codes, 0, sizeof(num_codes));
		for (uint32_t i = 0; i < num_syms; i++)
		{
			assert(pCode_sizes[i] <= FPNG_DECODER_TABLE_SIZE);
			num_codes[pCode_sizes[i]]++;
		}

		uint32_t next_code[17];
		next_code[0] = next_code[1] = 0;
		uint32_t total = 0;
		for (uint32_t i = 1; i <= 15; i++)
			next_code[i + 1] = (uint32_t)(total = ((total + ((uint32_t)num_codes[i])) << 1));

		if (total != 0x10000)
		{
			uint32_t j = 0;

			for (uint32_t i = 15; i != 0; i--)
				if ((j += num_codes[i]) > 1)
					return false;
			
			if (j != 1)
				return false;
		}

		uint32_t rev_codes[DEFL_MAX_HUFF_SYMBOLS];

		for (uint32_t i = 0; i < num_syms; i++)
			rev_codes[i] = next_code[pCode_sizes[i]]++;

		memset(pTable, 0, sizeof(uint32_t) * FPNG_DECODER_TABLE_SIZE);

		for (uint32_t i = 0; i < num_syms; i++)
		{
			const uint32_t code_size = pCode_sizes[i];
			if (!code_size)
				continue;

			uint32_t old_code = rev_codes[i], new_code = 0;
			for (uint32_t j = code_size; j != 0; j--)
			{
				new_code = (new_code << 1) | (old_code & 1);
				old_code >>= 1;
			}

			uint32_t j = 1 << code_size;

			while (new_code < FPNG_DECODER_TABLE_SIZE)
			{
				pTable[new_code] = i | (code_size << 9);
				new_code += j;
			}
		}

		return true;
	}

	static const uint16_t g_run_len3_to_4[259] = 
	{
		0,
		0, 0, 4, 0, 0, 8, 0, 0, 12, 0, 0, 16, 0, 0, 20, 0, 0, 24, 0, 0, 28, 0, 0,
		32, 0, 0, 36, 0, 0, 40, 0, 0, 44, 0, 0, 48, 0, 0, 52, 0, 0, 56, 0, 0,
		60, 0, 0, 64, 0, 0, 68, 0, 0, 72, 0, 0, 76, 0, 0, 80, 0, 0, 84, 0, 0,
		88, 0, 0, 92, 0, 0, 96, 0, 0, 100, 0, 0, 104, 0, 0, 108, 0, 0, 112, 0, 0,
		116, 0, 0, 120, 0, 0, 124, 0, 0, 128, 0, 0, 132, 0, 0, 136, 0, 0, 140, 0, 0,
		144, 0, 0, 148, 0, 0, 152, 0, 0, 156, 0, 0, 160, 0, 0, 164, 0, 0, 168, 0, 0,
		172, 0, 0, 176, 0, 0, 180, 0, 0, 184, 0, 0, 188, 0, 0, 192, 0, 0, 196, 0, 0,
		200, 0, 0, 204, 0, 0, 208, 0, 0, 212, 0, 0, 216, 0, 0, 220, 0, 0, 224, 0, 0,
		228, 0, 0, 232, 0, 0, 236, 0, 0, 240, 0, 0, 244, 0, 0, 248, 0, 0, 252, 0, 0,
		256, 0, 0, 260, 0, 0, 264, 0, 0, 268, 0, 0, 272, 0, 0, 276, 0, 0, 280, 0, 0,
		284, 0, 0, 288, 0, 0, 292, 0, 0, 296, 0, 0, 300, 0, 0, 304, 0, 0, 308, 0, 0,
		312, 0, 0, 316, 0, 0, 320, 0, 0, 324, 0, 0, 328, 0, 0, 332, 0, 0, 336, 0, 0,
		340, 0, 0, 
		344,
	};

	static const int s_length_extra[] = { 0,0,0,0, 0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, 0,    0,0 };
	static const int s_length_range[] = { 3,4,5,6, 7,8,9,10, 11,13,15,17, 19,23,27,31, 35,43,51,59, 67,83,99,115, 131,163,195,227, 258,    0,0 };

#define ENSURE_32BITS() do { \
	if (bit_buf_size < 32) { \
		if ((src_ofs + 4) > src_len) return false; \
		bit_buf |= ((uint64_t)READ_LE32(pSrc + src_ofs)) << bit_buf_size; \
		src_ofs += 4; bit_buf_size += 32; } \
	} while(0)

#define GET_BITS(b, ll) do { \
	uint32_t l = ll; assert(l && (l <= 32)); \
	b = (uint32_t)(bit_buf & g_bitmasks[l]); \
	bit_buf >>= l; \
	bit_buf_size -= l; \
	ENSURE_32BITS(); \
	} while(0)

#define SKIP_BITS(ll) do { \
	uint32_t l = ll; assert(l <= 32); \
	bit_buf >>= l; \
	bit_buf_size -= l; \
	ENSURE_32BITS(); \
	} while(0)

#define GET_BITS_NE(b, ll) do { \
	uint32_t l = ll; assert(l && (l <= 32) && (bit_buf_size >= l)); \
	b = (uint32_t)(bit_buf & g_bitmasks[l]); \
	bit_buf >>= l; \
	bit_buf_size -= l; \
	} while(0)

#define SKIP_BITS_NE(ll) do { \
	uint32_t l = ll; assert(l <= 32 && (bit_buf_size >= l)); \
	bit_buf >>= l; \
	bit_buf_size -= l; \
	} while(0)

	static bool prepare_dynamic_block(
		const uint8_t* pSrc, uint32_t src_len, uint32_t& src_ofs,
		uint32_t& bit_buf_size, uint64_t& bit_buf,
		uint32_t* pLit_table, uint32_t num_chans)
	{
		static const uint8_t s_bit_length_order[] = { 16, 17, 18, 0, 8,  7,  9, 6, 10,  5, 11, 4, 12,  3, 13, 2, 14,  1, 15 };

		uint32_t num_lit_codes, num_dist_codes, num_clen_codes;

		GET_BITS(num_lit_codes, 5);
		num_lit_codes += 257;

		GET_BITS(num_dist_codes, 5);
		num_dist_codes += 1;
		
		uint32_t total_codes = num_lit_codes + num_dist_codes;
		if (total_codes > (DEFL_MAX_HUFF_SYMBOLS_0 + DEFL_MAX_HUFF_SYMBOLS_1))
			return false;

		uint8_t code_sizes[DEFL_MAX_HUFF_SYMBOLS_0 + DEFL_MAX_HUFF_SYMBOLS_1];
		memset(code_sizes, 0, sizeof(code_sizes));

		GET_BITS(num_clen_codes, 4);
		num_clen_codes += 4;

		uint8_t clen_codesizes[DEFL_MAX_HUFF_SYMBOLS_2];
		memset(clen_codesizes, 0, sizeof(clen_codesizes));

		for (uint32_t i = 0; i < num_clen_codes; i++)
		{
			uint32_t len = 0;
			GET_BITS(len, 3);
			clen_codesizes[s_bit_length_order[i]] = (uint8_t)len;
		}

		uint32_t clen_table[FPNG_DECODER_TABLE_SIZE];
		if (!build_decoder_table(DEFL_MAX_HUFF_SYMBOLS_2, clen_codesizes, clen_table))
			return false;

		uint32_t min_code_size = 15;

		for (uint32_t cur_code = 0; cur_code < total_codes; )
		{
			uint32_t sym = clen_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
			uint32_t sym_len = sym >> 9;
			if (!sym_len)
				return false;
			SKIP_BITS(sym_len);
			sym &= 511;
						
			if (sym <= 15)
			{
				// Can't be a fpng Huffman table
				if (sym > FPNG_DECODER_TABLE_BITS)
					return false;

				if (sym)
					min_code_size = minimum(min_code_size, sym);

				code_sizes[cur_code++] = (uint8_t)sym;
				continue;
			}

			uint32_t rep_len = 0, rep_code_size = 0;

			switch (sym)
			{
			case 16:
			{
				GET_BITS(rep_len, 2);
				rep_len += 3;
				if (!cur_code)
					return false;
				rep_code_size = code_sizes[cur_code - 1];
				break;
			}
			case 17:
			{
				GET_BITS(rep_len, 3);
				rep_len += 3;
				rep_code_size = 0;
				break;
			}
			case 18:
			{
				GET_BITS(rep_len, 7);
				rep_len += 11;
				rep_code_size = 0;
				break;
			}
			}

			if ((cur_code + rep_len) > total_codes)
				return false;

			for (; rep_len; rep_len--)
				code_sizes[cur_code++] = (uint8_t)rep_code_size;
		}

		uint8_t lit_codesizes[DEFL_MAX_HUFF_SYMBOLS_0];

		memcpy(lit_codesizes, code_sizes, num_lit_codes);
		memset(lit_codesizes + num_lit_codes, 0, DEFL_MAX_HUFF_SYMBOLS_0 - num_lit_codes);

		uint32_t total_valid_distcodes = 0;
		for (uint32_t i = 0; i < num_dist_codes; i++)
			total_valid_distcodes += (code_sizes[num_lit_codes + i] == 1);
		
		// 1 or 2 because the first version of FPNG only issued 1 valid distance code, but that upset wuffs. So we let 1 or 2 through.
		if ((total_valid_distcodes < 1) || (total_valid_distcodes > 2))
			return false;

		if (code_sizes[num_lit_codes + (num_chans - 1)] != 1)
			return false;

		if (total_valid_distcodes == 2)
		{
			// If there are two valid distance codes, make sure the first is 1 bit.
			if (code_sizes[num_lit_codes + num_chans] != 1)
				return false;
		}
						
		if (!build_decoder_table(num_lit_codes, lit_codesizes, pLit_table))
			return false;

		// Add next symbol to decoder table, when it fits
		for (uint32_t i = 0; i < FPNG_DECODER_TABLE_SIZE; i++)
		{
			uint32_t sym = pLit_table[i] & 511;
			if (sym >= 256)
				continue;

			uint32_t sym_bits = (pLit_table[i] >> 9) & 15;
			if (!sym_bits)
				continue;
			assert(sym_bits <= FPNG_DECODER_TABLE_BITS);

			uint32_t bits_left = FPNG_DECODER_TABLE_BITS - sym_bits;
			if (bits_left < min_code_size)
				continue;

			uint32_t next_bits = i >> sym_bits;
			uint32_t next_sym = pLit_table[next_bits] & 511;
			uint32_t next_sym_bits = (pLit_table[next_bits] >> 9) & 15;
			if ((!next_sym_bits) || (bits_left < next_sym_bits))
				continue;

			pLit_table[i] |= (next_sym << 16) | (next_sym_bits << (16 + 9));
		}

		return true;
	}
		
	static bool fpng_pixel_zlib_raw_decompress(
		const uint8_t* pSrc, uint32_t src_len, uint32_t zlib_len,
		uint8_t* pDst, uint32_t w, uint32_t h,
		uint32_t src_chans, uint32_t dst_chans)
	{
		assert((src_chans == 3) || (src_chans == 4));
		assert((dst_chans == 3) || (dst_chans == 4));
		
		const uint32_t src_bpl = w * src_chans;
		const uint32_t dst_bpl = w * dst_chans;
		const uint32_t dst_len = dst_bpl * h;

		uint32_t src_ofs = 2;
		uint32_t dst_ofs = 0;
		uint32_t raster_ofs = 0;
		uint32_t comp_ofs = 0;

		for (; ; )
		{
			if ((src_ofs + 1) > src_len)
				return false;

			const bool bfinal = (pSrc[src_ofs] & 1) != 0;
			const uint32_t btype = (pSrc[src_ofs] >> 1) & 3;
			if (btype != 0)
				return false;

			src_ofs++;

			if ((src_ofs + 4) > src_len)
				return false;
			uint32_t len = pSrc[src_ofs + 0] | (pSrc[src_ofs + 1] << 8);
			uint32_t nlen = pSrc[src_ofs + 2] | (pSrc[src_ofs + 3] << 8);
			src_ofs += 4;

			if (len != (~nlen & 0xFFFF))
				return false;

			if ((src_ofs + len) > src_len)
				return false;

			// Raw blocks are a relatively uncommon case so this isn't well optimized.
			// Supports 3->4 and 4->3 byte/pixel conversion.
			for (uint32_t i = 0; i < len; i++)
			{
				uint32_t c = pSrc[src_ofs + i];

				if (!raster_ofs)
				{
					// Check filter type
					if (c != 0)
						return false;
					
					assert(!comp_ofs);
				}
				else
				{
					if (comp_ofs < dst_chans)
					{
						if (dst_ofs == dst_len)
							return false;

						pDst[dst_ofs++] = (uint8_t)c;
					}
					
					if (++comp_ofs == src_chans)
					{
						if (dst_chans > src_chans)
						{
							if (dst_ofs == dst_len)
								return false;

							pDst[dst_ofs++] = (uint8_t)0xFF;
						}

						comp_ofs = 0;
					}
				}

				if (++raster_ofs == (src_bpl + 1))
				{
					assert(!comp_ofs);
					raster_ofs = 0;
				}
			}

			src_ofs += len;

			if (bfinal)
				break;
		}

		if (comp_ofs != 0)
			return false;

		// Check for zlib adler32
		if ((src_ofs + 4) != zlib_len)
			return false;

		return (dst_ofs == dst_len);
	}
	
	template<uint32_t dst_comps>
	static bool fpng_pixel_zlib_decompress_3(
		const uint8_t* pSrc, uint32_t src_len, uint32_t zlib_len,
		uint8_t* pDst, uint32_t w, uint32_t h)
	{
		assert(src_len >= (zlib_len + 4));

		const uint32_t dst_bpl = w * dst_comps;
		//const uint32_t dst_len = dst_bpl * h;

		if (zlib_len < 7)
			return false;

		// check zlib header
		if ((pSrc[0] != 0x78) || (pSrc[1] != 0x01))
			return false;

		uint32_t src_ofs = 2;
		
		if ((pSrc[src_ofs] & 6) == 0)
			return fpng_pixel_zlib_raw_decompress(pSrc, src_len, zlib_len, pDst, w, h, 3, dst_comps);
		
		if ((src_ofs + 4) > src_len)
			return false;
		uint64_t bit_buf = READ_LE32(pSrc + src_ofs);
		src_ofs += 4;

		uint32_t bit_buf_size = 32;

		uint32_t bfinal, btype;
		GET_BITS(bfinal, 1);
		GET_BITS(btype, 2);

		// Must be the final block or it's not valid, and type=2 (dynamic)
		if ((bfinal != 1) || (btype != 2))
			return false;
		
		uint32_t lit_table[FPNG_DECODER_TABLE_SIZE];
		if (!prepare_dynamic_block(pSrc, src_len, src_ofs, bit_buf_size, bit_buf, lit_table, 3))
			return false;

		const uint8_t* pPrev_scanline = nullptr;
		uint8_t* pCur_scanline = pDst;

		for (uint32_t y = 0; y < h; y++)
		{
			// At start of PNG scanline, so read the filter literal
			assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
			uint32_t filter = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
			uint32_t filter_len = (filter >> 9) & 15;
			if (!filter_len)
				return false;
			SKIP_BITS(filter_len);
			filter &= 511;

			uint32_t expected_filter = (y ? 2 : 0);
			if (filter != expected_filter)
				return false;

			uint32_t x_ofs = 0;
			uint8_t prev_delta_r = 0, prev_delta_g = 0, prev_delta_b = 0;
			do
			{
				assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
				uint32_t lit0_tab = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
				
				uint32_t lit0 = lit0_tab;
				uint32_t lit0_len = (lit0_tab >> 9) & 15;
				if (!lit0_len)
					return false;
				SKIP_BITS(lit0_len);

				if (lit0 & 256)
				{
					lit0 &= 511;

					// Can't be EOB - we still have more pixels to decompress.
					if (lit0 == 256)
						return false;

					// Must be an RLE match against the previous pixel.
					uint32_t run_len = s_length_range[lit0 - 257];
					if (lit0 >= 265)
					{
						uint32_t e;
						GET_BITS_NE(e, s_length_extra[lit0 - 257]);

						run_len += e;
					}
					
					// Skip match distance - it's always the same (3)
					SKIP_BITS_NE(1);

					// Matches must always be a multiple of 3/4 bytes
					assert((run_len % 3) == 0);
																				
					if (dst_comps == 4)
					{
						const uint32_t x_ofs_end = x_ofs + g_run_len3_to_4[run_len];
						
						// Check for valid run lengths
						if (x_ofs == x_ofs_end)
							return false;

						// Matches cannot cross scanlines.
						if (x_ofs_end > dst_bpl)
							return false;

						if (pPrev_scanline)
						{
							if ((prev_delta_r | prev_delta_g | prev_delta_b) == 0)
							{
								memcpy(pCur_scanline + x_ofs, pPrev_scanline + x_ofs, x_ofs_end - x_ofs);
								x_ofs = x_ofs_end;
							}
							else
							{
								do
								{
									pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + prev_delta_r);
									pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + prev_delta_g);
									pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + prev_delta_b);
									pCur_scanline[x_ofs + 3] = 0xFF;
									x_ofs += 4;
								} while (x_ofs < x_ofs_end);
							}
						}
						else
						{
							do
							{
								pCur_scanline[x_ofs] = prev_delta_r;
								pCur_scanline[x_ofs + 1] = prev_delta_g;
								pCur_scanline[x_ofs + 2] = prev_delta_b;
								pCur_scanline[x_ofs + 3] = 0xFF;
								x_ofs += 4;
							} while (x_ofs < x_ofs_end);
						}
					}
					else
					{
						// Check for valid run lengths
						if (!g_run_len3_to_4[run_len])
							return false;

						const uint32_t x_ofs_end = x_ofs + run_len;

						// Matches cannot cross scanlines.
						if (x_ofs_end > dst_bpl)
							return false;

						if (pPrev_scanline)
						{
							if ((prev_delta_r | prev_delta_g | prev_delta_b) == 0)
							{
								memcpy(pCur_scanline + x_ofs, pPrev_scanline + x_ofs, run_len);
								x_ofs = x_ofs_end;
							}
							else
							{
								do
								{
									pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + prev_delta_r);
									pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + prev_delta_g);
									pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + prev_delta_b);
									x_ofs += 3;
								} while (x_ofs < x_ofs_end);
							}
						}
						else
						{
							do
							{
								pCur_scanline[x_ofs] = prev_delta_r;
								pCur_scanline[x_ofs + 1] = prev_delta_g;
								pCur_scanline[x_ofs + 2] = prev_delta_b;
								x_ofs += 3;
							} while (x_ofs < x_ofs_end);
						}
					}
				}
				else
				{
					uint32_t lit1, lit2;

					uint32_t lit1_spec_len = (lit0_tab >> (16 + 9));
					uint32_t lit2_len;
					if (lit1_spec_len)
					{
						lit1 = (lit0_tab >> 16) & 511;
						SKIP_BITS_NE(lit1_spec_len);

						assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
						lit2 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
						lit2_len = (lit2 >> 9) & 15;
						if (!lit2_len)
							return false;
					}
					else
					{
						assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
						lit1 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
						uint32_t lit1_len = (lit1 >> 9) & 15;
						if (!lit1_len)
							return false;
						SKIP_BITS_NE(lit1_len);

						lit2_len = (lit1 >> (16 + 9));
						if (lit2_len)
							lit2 = lit1 >> 16;
						else
						{
							assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
							lit2 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
							lit2_len = (lit2 >> 9) & 15;
							if (!lit2_len)
								return false;
						}
					}

					SKIP_BITS(lit2_len);
					
					// Check for matches
					if ((lit1 | lit2) & 256)
						return false;

					if (dst_comps == 4)
					{
						if (pPrev_scanline)
						{
							pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + lit0);
							pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + lit1);
							pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + lit2);
							pCur_scanline[x_ofs + 3] = 0xFF;
						}
						else
						{
							pCur_scanline[x_ofs] = (uint8_t)lit0;
							pCur_scanline[x_ofs + 1] = (uint8_t)lit1;
							pCur_scanline[x_ofs + 2] = (uint8_t)lit2;
							pCur_scanline[x_ofs + 3] = 0xFF;
						}
						x_ofs += 4;
					}
					else
					{
						if (pPrev_scanline)
						{
							pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + lit0);
							pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + lit1);
							pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + lit2);
						}
						else
						{
							pCur_scanline[x_ofs] = (uint8_t)lit0;
							pCur_scanline[x_ofs + 1] = (uint8_t)lit1;
							pCur_scanline[x_ofs + 2] = (uint8_t)lit2;
						}
						x_ofs += 3;
					}

					prev_delta_r = (uint8_t)lit0;
					prev_delta_g = (uint8_t)lit1;
					prev_delta_b = (uint8_t)lit2;
										
					// See if we can decode one more pixel.
					uint32_t spec_next_len0_len = lit2 >> (16 + 9);
					if ((spec_next_len0_len) && (x_ofs < dst_bpl))
					{
						lit0 = (lit2 >> 16) & 511;
						if (lit0 < 256)
						{
							SKIP_BITS_NE(spec_next_len0_len);

							assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
							lit1 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
							uint32_t lit1_len = (lit1 >> 9) & 15;
							if (!lit1_len)
								return false;
							SKIP_BITS(lit1_len);

							lit2_len = (lit1 >> (16 + 9));
							if (lit2_len)
								lit2 = lit1 >> 16;
							else
							{
								assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
								lit2 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
								lit2_len = (lit2 >> 9) & 15;
								if (!lit2_len)
									return false;
							}

							SKIP_BITS_NE(lit2_len);

							// Check for matches
							if ((lit1 | lit2) & 256)
								return false;
					
							if (dst_comps == 4)
							{
								if (pPrev_scanline)
								{
									pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + lit0);
									pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + lit1);
									pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + lit2);
									pCur_scanline[x_ofs + 3] = 0xFF;
								}
								else
								{
									pCur_scanline[x_ofs] = (uint8_t)lit0;
									pCur_scanline[x_ofs + 1] = (uint8_t)lit1;
									pCur_scanline[x_ofs + 2] = (uint8_t)lit2;
									pCur_scanline[x_ofs + 3] = 0xFF;
								}
								x_ofs += 4;
							}
							else
							{
								if (pPrev_scanline)
								{
									pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + lit0);
									pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + lit1);
									pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + lit2);
								}
								else
								{
									pCur_scanline[x_ofs] = (uint8_t)lit0;
									pCur_scanline[x_ofs + 1] = (uint8_t)lit1;
									pCur_scanline[x_ofs + 2] = (uint8_t)lit2;
								}
								x_ofs += 3;
							}

							prev_delta_r = (uint8_t)lit0;
							prev_delta_g = (uint8_t)lit1;
							prev_delta_b = (uint8_t)lit2;
																				
						} // if (lit0 < 256)

					} // if ((spec_next_len0_len) && (x_ofs < bpl))
				}

			} while (x_ofs < dst_bpl);

			pPrev_scanline = pCur_scanline;
			pCur_scanline += dst_bpl;

		} // y

		// The last symbol should be EOB
		assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
		uint32_t lit0 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
		uint32_t lit0_len = (lit0 >> 9) & 15;
		if (!lit0_len)
			return false;
		lit0 &= 511;
		if (lit0 != 256)
			return false;

		bit_buf_size -= lit0_len;
		bit_buf >>= lit0_len;

		uint32_t align_bits = bit_buf_size & 7;
		bit_buf_size -= align_bits;
		bit_buf >>= align_bits;

		if (src_ofs < (bit_buf_size >> 3))
			return false;
		src_ofs -= (bit_buf_size >> 3);

		// We should be at the very end, because the bit buf reads ahead 32-bits (which contains the zlib adler32).
		if ((src_ofs + 4) != zlib_len)
			return false;

		return true;
	}

	template<uint32_t dst_comps>
	static bool fpng_pixel_zlib_decompress_4(
		const uint8_t* pSrc, uint32_t src_len, uint32_t zlib_len,
		uint8_t* pDst, uint32_t w, uint32_t h)
	{
		assert(src_len >= (zlib_len + 4));

		const uint32_t dst_bpl = w * dst_comps;
		//const uint32_t dst_len = dst_bpl * h;

		if (zlib_len < 7)
			return false;

		// check zlib header
		if ((pSrc[0] != 0x78) || (pSrc[1] != 0x01))
			return false;

		uint32_t src_ofs = 2;

		if ((pSrc[src_ofs] & 6) == 0)
			return fpng_pixel_zlib_raw_decompress(pSrc, src_len, zlib_len, pDst, w, h, 4, dst_comps);

		if ((src_ofs + 4) > src_len)
			return false;
		uint64_t bit_buf = READ_LE32(pSrc + src_ofs);
		src_ofs += 4;

		uint32_t bit_buf_size = 32;

		uint32_t bfinal, btype;
		GET_BITS(bfinal, 1);
		GET_BITS(btype, 2);

		// Must be the final block or it's not valid, and type=2 (dynamic)
		if ((bfinal != 1) || (btype != 2))
			return false;

		uint32_t lit_table[FPNG_DECODER_TABLE_SIZE];
		if (!prepare_dynamic_block(pSrc, src_len, src_ofs, bit_buf_size, bit_buf, lit_table, 4))
			return false;

		const uint8_t* pPrev_scanline = nullptr;
		uint8_t* pCur_scanline = pDst;

		for (uint32_t y = 0; y < h; y++)
		{
			// At start of PNG scanline, so read the filter literal
			assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
			uint32_t filter = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
			uint32_t filter_len = (filter >> 9) & 15;
			if (!filter_len)
				return false;
			SKIP_BITS(filter_len);
			filter &= 511;

			uint32_t expected_filter = (y ? 2 : 0);
			if (filter != expected_filter)
				return false;

			uint32_t x_ofs = 0;
			uint8_t prev_delta_r = 0, prev_delta_g = 0, prev_delta_b = 0, prev_delta_a = 0;
			do
			{
				assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
				uint32_t lit0_tab = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];

				uint32_t lit0 = lit0_tab;
				uint32_t lit0_len = (lit0_tab >> 9) & 15;
				if (!lit0_len)
					return false;
				SKIP_BITS(lit0_len);

				if (lit0 & 256)
				{
					lit0 &= 511;

					// Can't be EOB - we still have more pixels to decompress.
					if (lit0 == 256)
						return false;

					// Must be an RLE match against the previous pixel.
					uint32_t run_len = s_length_range[lit0 - 257];
					if (lit0 >= 265)
					{
						uint32_t e;
						GET_BITS_NE(e, s_length_extra[lit0 - 257]);

						run_len += e;
					}

					// Skip match distance - it's always the same (4)
					SKIP_BITS_NE(1);

					// Matches must always be a multiple of 3/4 bytes
					if (run_len & 3)
						return false;
										
					if (dst_comps == 3)
					{
						const uint32_t run_len3 = (run_len >> 2) * 3;
						const uint32_t x_ofs_end = x_ofs + run_len3;

						// Matches cannot cross scanlines.
						if (x_ofs_end > dst_bpl)
							return false;

						if (pPrev_scanline)
						{
							if ((prev_delta_r | prev_delta_g | prev_delta_b | prev_delta_a) == 0)
							{
								memcpy(pCur_scanline + x_ofs, pPrev_scanline + x_ofs, run_len3);
								x_ofs = x_ofs_end;
							}
							else
							{
								do
								{
									pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + prev_delta_r);
									pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + prev_delta_g);
									pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + prev_delta_b);
									x_ofs += 3;
								} while (x_ofs < x_ofs_end);
							}
						}
						else
						{
							do
							{
								pCur_scanline[x_ofs] = prev_delta_r;
								pCur_scanline[x_ofs + 1] = prev_delta_g;
								pCur_scanline[x_ofs + 2] = prev_delta_b;
								x_ofs += 3;
							} while (x_ofs < x_ofs_end);
						}
					}
					else
					{
						const uint32_t x_ofs_end = x_ofs + run_len;

						// Matches cannot cross scanlines.
						if (x_ofs_end > dst_bpl)
							return false;

						if (pPrev_scanline)
						{
							if ((prev_delta_r | prev_delta_g | prev_delta_b | prev_delta_a) == 0)
							{
								memcpy(pCur_scanline + x_ofs, pPrev_scanline + x_ofs, run_len);
								x_ofs = x_ofs_end;
							}
							else
							{
								do
								{
									pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + prev_delta_r);
									pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + prev_delta_g);
									pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + prev_delta_b);
									pCur_scanline[x_ofs + 3] = (uint8_t)(pPrev_scanline[x_ofs + 3] + prev_delta_a);
									x_ofs += 4;
								} while (x_ofs < x_ofs_end);
							}
						}
						else
						{
							do
							{
								pCur_scanline[x_ofs] = prev_delta_r;
								pCur_scanline[x_ofs + 1] = prev_delta_g;
								pCur_scanline[x_ofs + 2] = prev_delta_b;
								pCur_scanline[x_ofs + 3] = prev_delta_a;
								x_ofs += 4;
							} while (x_ofs < x_ofs_end);
						}
					}
				}
				else
				{
					uint32_t lit1, lit2;

					uint32_t lit1_spec_len = (lit0_tab >> (16 + 9));
					uint32_t lit2_len;
					if (lit1_spec_len)
					{
						lit1 = (lit0_tab >> 16) & 511;
						SKIP_BITS_NE(lit1_spec_len);

						assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
						lit2 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
						lit2_len = (lit2 >> 9) & 15;
						if (!lit2_len)
							return false;
					}
					else
					{
						assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
						lit1 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
						uint32_t lit1_len = (lit1 >> 9) & 15;
						if (!lit1_len)
							return false;
						SKIP_BITS_NE(lit1_len);

						lit2_len = (lit1 >> (16 + 9));
						if (lit2_len)
							lit2 = lit1 >> 16;
						else
						{
							assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
							lit2 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
							lit2_len = (lit2 >> 9) & 15;
							if (!lit2_len)
								return false;
						}
					}

					uint32_t lit3;
					uint32_t lit3_len = lit2 >> (16 + 9);
					
					if (lit3_len)
					{
						lit3 = (lit2 >> 16);
						SKIP_BITS(lit2_len + lit3_len);
					}
					else
					{
						SKIP_BITS(lit2_len);

						assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
						lit3 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
						lit3_len = (lit3 >> 9) & 15;
						if (!lit3_len)
							return false;

						SKIP_BITS_NE(lit3_len);
					}
										
					// Check for matches
					if ((lit1 | lit2 | lit3) & 256)
						return false;

					if (dst_comps == 3)
					{
						if (pPrev_scanline)
						{
							pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + lit0);
							pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + lit1);
							pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + lit2);
						}
						else
						{
							pCur_scanline[x_ofs] = (uint8_t)lit0;
							pCur_scanline[x_ofs + 1] = (uint8_t)lit1;
							pCur_scanline[x_ofs + 2] = (uint8_t)lit2;
						}

						x_ofs += 3;
					}
					else
					{
						if (pPrev_scanline)
						{
							pCur_scanline[x_ofs] = (uint8_t)(pPrev_scanline[x_ofs] + lit0);
							pCur_scanline[x_ofs + 1] = (uint8_t)(pPrev_scanline[x_ofs + 1] + lit1);
							pCur_scanline[x_ofs + 2] = (uint8_t)(pPrev_scanline[x_ofs + 2] + lit2);
							pCur_scanline[x_ofs + 3] = (uint8_t)(pPrev_scanline[x_ofs + 3] + lit3);
						}
						else
						{
							pCur_scanline[x_ofs] = (uint8_t)lit0;
							pCur_scanline[x_ofs + 1] = (uint8_t)lit1;
							pCur_scanline[x_ofs + 2] = (uint8_t)lit2;
							pCur_scanline[x_ofs + 3] = (uint8_t)lit3;
						}
						
						x_ofs += 4;
					}

					prev_delta_r = (uint8_t)lit0;
					prev_delta_g = (uint8_t)lit1;
					prev_delta_b = (uint8_t)lit2;
					prev_delta_a = (uint8_t)lit3;
				}

			} while (x_ofs < dst_bpl);

			pPrev_scanline = pCur_scanline;
			pCur_scanline += dst_bpl;
		} // y

		// The last symbol should be EOB
		assert(bit_buf_size >= FPNG_DECODER_TABLE_BITS);
		uint32_t lit0 = lit_table[bit_buf & (FPNG_DECODER_TABLE_SIZE - 1)];
		uint32_t lit0_len = (lit0 >> 9) & 15;
		if (!lit0_len)
			return false;
		lit0 &= 511;
		if (lit0 != 256)
			return false;

		bit_buf_size -= lit0_len;
		bit_buf >>= lit0_len;

		uint32_t align_bits = bit_buf_size & 7;
		bit_buf_size -= align_bits;
		bit_buf >>= align_bits;

		if (src_ofs < (bit_buf_size >> 3))
			return false;
		src_ofs -= (bit_buf_size >> 3);

		// We should be at the very end, because the bit buf reads ahead 32-bits (which contains the zlib adler32).
		if ((src_ofs + 4) != zlib_len)
			return false;

		return true;
	}

#pragma pack(push)
#pragma pack(1)
	struct png_chunk_prefix
	{
		uint32_t m_length;
		uint8_t m_type[4];
	};
	struct png_ihdr
	{
		png_chunk_prefix m_prefix;
		uint32_t m_width;
		uint32_t m_height;
		uint8_t m_bitdepth;
		uint8_t m_color_type;
		uint8_t m_comp_method;
		uint8_t m_filter_method;
		uint8_t m_interlace_method;
		uint32_t m_crc32;
	};
	const uint32_t IHDR_EXPECTED_LENGTH = 13;
	struct png_iend
	{
		png_chunk_prefix m_prefix;
		uint32_t m_crc32;
	};
#pragma pack(pop)

	static int fpng_get_info_internal(const void* pImage, uint32_t image_size, uint32_t& width, uint32_t& height, uint32_t& channels_in_file, uint32_t &idat_ofs, uint32_t &idat_len)
	{
		static const uint8_t s_png_sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };

		if (!endian_check())
		{
			assert(0);
			return false;
		}
				
		width = 0;
		height = 0;
		channels_in_file = 0;
		idat_ofs = 0, idat_len = 0;
				
		// Ensure the file has at least a minimum possible size
		if (image_size < (sizeof(s_png_sig) + sizeof(png_ihdr) + sizeof(png_chunk_prefix) + 1 + sizeof(uint32_t) + sizeof(png_iend)))
			return FPNG_DECODE_FAILED_NOT_PNG;

		if (memcmp(pImage, s_png_sig, 8) != 0)
			return FPNG_DECODE_FAILED_NOT_PNG;

		const uint8_t* pImage_u8 = static_cast<const uint8_t*>(pImage) + 8;

		const png_ihdr& ihdr = *reinterpret_cast<const png_ihdr*>(pImage_u8);
		pImage_u8 += sizeof(png_ihdr);

		if (READ_BE32(&ihdr.m_prefix.m_length) != IHDR_EXPECTED_LENGTH)
			return FPNG_DECODE_FAILED_NOT_PNG;

		if (fpng_crc32(ihdr.m_prefix.m_type, 4 + IHDR_EXPECTED_LENGTH, FPNG_CRC32_INIT) != READ_BE32(&ihdr.m_crc32))
			return FPNG_DECODE_FAILED_HEADER_CRC32;

		width = READ_BE32(&ihdr.m_width);
		height = READ_BE32(&ihdr.m_height);
				
		if (!width || !height || (width > FPNG_MAX_SUPPORTED_DIM) || (height > FPNG_MAX_SUPPORTED_DIM))
			return FPNG_DECODE_FAILED_INVALID_DIMENSIONS;

		uint64_t total_pixels = (uint64_t)width * height;
		if (total_pixels > (1 << 30))
			return FPNG_DECODE_FAILED_INVALID_DIMENSIONS;

		if ((ihdr.m_comp_method) || (ihdr.m_filter_method) || (ihdr.m_interlace_method) || (ihdr.m_bitdepth != 8))
			return FPNG_DECODE_NOT_FPNG;

		if (ihdr.m_color_type == 2)
			channels_in_file = 3;
		else if (ihdr.m_color_type == 6)
			channels_in_file = 4;

		if (!channels_in_file)
			return FPNG_DECODE_NOT_FPNG;

		// Scan all the chunks. Look for one IDAT, IEND, and our custom fdEC chunk that indicates the file was compressed by us. Skip any ancillary chunks.
		bool found_fdec_chunk = false;
		
		for (; ; )
		{
			const size_t src_ofs = pImage_u8 - static_cast<const uint8_t*>(pImage);
			if (src_ofs >= image_size)
				return FPNG_DECODE_FAILED_CHUNK_PARSING;

			const uint32_t bytes_remaining = image_size - (uint32_t)src_ofs;
			if (bytes_remaining < sizeof(uint32_t) * 3)
				return FPNG_DECODE_FAILED_CHUNK_PARSING;

			const png_chunk_prefix* pChunk = reinterpret_cast<const png_chunk_prefix*>(pImage_u8);

			const uint32_t chunk_len = READ_BE32(&pChunk->m_length);
			if ((src_ofs + sizeof(uint32_t) + chunk_len + sizeof(uint32_t)) > image_size)
				return FPNG_DECODE_FAILED_CHUNK_PARSING;

			for (uint32_t i = 0; i < 4; i++)
			{
				const uint8_t c = pChunk->m_type[i];
				const bool is_upper = (c >= 65) && (c <= 90), is_lower = (c >= 97) && (c <= 122);
				if ((!is_upper) && (!is_lower))
					return FPNG_DECODE_FAILED_CHUNK_PARSING;
			}

			const uint32_t expected_crc32 = READ_BE32(pImage_u8 + sizeof(uint32_t) * 2 + chunk_len);

			char chunk_type[5] = { (char)pChunk->m_type[0], (char)pChunk->m_type[1], (char)pChunk->m_type[2], (char)pChunk->m_type[3], 0 };
			const bool is_idat = strcmp(chunk_type, "IDAT") == 0;

#if !FPNG_DISABLE_DECODE_CRC32_CHECKS
			if (!is_idat)
			{
				uint32_t actual_crc32 = fpng_crc32(pImage_u8 + sizeof(uint32_t), sizeof(uint32_t) + chunk_len, FPNG_CRC32_INIT);
				if (actual_crc32 != expected_crc32)
					return FPNG_DECODE_FAILED_HEADER_CRC32;
			}
#endif

			const uint8_t* pChunk_data = pImage_u8 + sizeof(uint32_t) * 2;

			if (strcmp(chunk_type, "IEND") == 0)
				break;
			else if (is_idat)
			{
				// If there were multiple IDAT's, or we didn't find the fdEC chunk, then it's not FPNG.
				if ((idat_ofs) || (!found_fdec_chunk))
					return FPNG_DECODE_NOT_FPNG;

				idat_ofs = (uint32_t)src_ofs;
				idat_len = chunk_len;

				// Sanity check the IDAT chunk length
				if (idat_len < 7)
					return FPNG_DECODE_FAILED_INVALID_IDAT;
			}
			else if (strcmp(chunk_type, "fdEC") == 0)
			{
				if (found_fdec_chunk)
					return FPNG_DECODE_NOT_FPNG;

				// We've got our fdEC chunk. Now make sure it's big enough and check its contents.
				if (chunk_len != 5)
					return FPNG_DECODE_NOT_FPNG;

				// Check fdEC chunk sig
				if ((pChunk_data[0] != 82) || (pChunk_data[1] != 36) || (pChunk_data[2] != 147) || (pChunk_data[3] != 227))
					return FPNG_DECODE_NOT_FPNG;

				// Check fdEC version
				if (pChunk_data[4] != FPNG_FDEC_VERSION)
					return FPNG_DECODE_NOT_FPNG;

				found_fdec_chunk = true;
			}
			else
			{
				// Bail if it's a critical chunk - can't be FPNG
				if ((chunk_type[0] & 32) == 0)
					return FPNG_DECODE_NOT_FPNG;

				// ancillary chunk - skip it
			}

			pImage_u8 += sizeof(png_chunk_prefix) + chunk_len + sizeof(uint32_t);
		}

		if ((!found_fdec_chunk) || (!idat_ofs))
			return FPNG_DECODE_NOT_FPNG;
		
		return FPNG_DECODE_SUCCESS;
	}

	int fpng_get_info(const void* pImage, uint32_t image_size, uint32_t& width, uint32_t& height, uint32_t& channels_in_file)
	{
		uint32_t idat_ofs = 0, idat_len = 0;
		return fpng_get_info_internal(pImage, image_size, width, height, channels_in_file, idat_ofs, idat_len);
	}

	int fpng_decode_memory(const void *pImage, uint32_t image_size, std::vector<uint8_t> &out, uint32_t& width, uint32_t& height, uint32_t &channels_in_file, uint32_t desired_channels)
	{
		out.resize(0);
		width = 0;
		height = 0;
		channels_in_file = 0;

		if ((!pImage) || (!image_size) || ((desired_channels != 3) && (desired_channels != 4)))
		{
			assert(0);
			return FPNG_DECODE_INVALID_ARG;
		}

		uint32_t idat_ofs = 0, idat_len = 0;
		int status = fpng_get_info_internal(pImage, image_size, width, height, channels_in_file, idat_ofs, idat_len);
		if (status)
			return status;
				
		const uint64_t mem_needed = (uint64_t)width * height * desired_channels;
		if (mem_needed > UINT32_MAX)
			return FPNG_DECODE_FAILED_DIMENSIONS_TOO_LARGE;

		// On 32-bit systems do a quick sanity check before we try to resize the output buffer.
		if ((sizeof(size_t) == sizeof(uint32_t)) && (mem_needed >= 0x80000000))
			return FPNG_DECODE_FAILED_DIMENSIONS_TOO_LARGE;

		out.resize(mem_needed);
		
		const uint8_t* pIDAT_data = static_cast<const uint8_t*>(pImage) + idat_ofs + sizeof(uint32_t) * 2;
		const uint32_t src_len = image_size - (idat_ofs + sizeof(uint32_t) * 2);

		bool decomp_status;
		if (desired_channels == 3)
		{
			if (channels_in_file == 3)
				decomp_status = fpng_pixel_zlib_decompress_3<3>(pIDAT_data, src_len, idat_len, out.data(), width, height);
			else
				decomp_status = fpng_pixel_zlib_decompress_4<3>(pIDAT_data, src_len, idat_len, out.data(), width, height);
		}
		else
		{
			if (channels_in_file == 3)
				decomp_status = fpng_pixel_zlib_decompress_3<4>(pIDAT_data, src_len, idat_len, out.data(), width, height);
			else
				decomp_status = fpng_pixel_zlib_decompress_4<4>(pIDAT_data, src_len, idat_len, out.data(), width, height);
		}
		if (!decomp_status)
		{
			// Something went wrong. Either the file data was corrupted, or it doesn't conform to one of our zlib/Deflate constraints.
			// The conservative thing to do is indicate it wasn't written by us, and let the general purpose PNG decoder handle it.
			return FPNG_DECODE_NOT_FPNG;
		}

		return FPNG_DECODE_SUCCESS;
	}

#ifndef FPNG_NO_STDIO
	int fpng_decode_file(const char* pFilename, std::vector<uint8_t>& out, uint32_t& width, uint32_t& height, uint32_t& channels_in_file, uint32_t desired_channels)
	{
		FILE* pFile = nullptr;

#ifdef _MSC_VER
		fopen_s(&pFile, pFilename, "rb");
#else
		pFile = fopen(pFilename, "rb");
#endif

		if (!pFile)
			return FPNG_DECODE_FILE_OPEN_FAILED;

		if (fseek(pFile, 0, SEEK_END) != 0)
		{
			fclose(pFile);
			return FPNG_DECODE_FILE_SEEK_FAILED;
		}

#ifdef _WIN32
		int64_t filesize = _ftelli64(pFile);
#else
		int64_t filesize = ftello(pFile);
#endif

		if (fseek(pFile, 0, SEEK_SET) != 0)
		{
			fclose(pFile);
			return FPNG_DECODE_FILE_SEEK_FAILED;
		}

		if ( (filesize < 0) || (filesize > UINT32_MAX) || ( (sizeof(size_t) == sizeof(uint32_t)) && (filesize > 0x70000000) ) )
		{
			fclose(pFile);
			return FPNG_DECODE_FILE_TOO_LARGE;
		}

		std::vector<uint8_t> buf((size_t)filesize);
		if (fread(buf.data(), 1, buf.size(), pFile) != buf.size())
		{
			fclose(pFile);
			return FPNG_DECODE_FILE_READ_FAILED;
		}

		fclose(pFile);

		return fpng_decode_memory(buf.data(), (uint32_t)buf.size(), out, width, height, channels_in_file, desired_channels);
	}
#endif

} // namespace fpng

/*
	This is free and unencumbered software released into the public domain.

	Anyone is free to copy, modify, publish, use, compile, sell, or
	distribute this software, either in source code form or as a compiled
	binary, for any purpose, commercial or non-commercial, and by any
	means.

	In jurisdictions that recognize copyright laws, the author or authors
	of this software dedicate any and all copyright interest in the
	software to the public domain. We make this dedication for the benefit
	of the public at large and to the detriment of our heirs and
	successors. We intend this dedication to be an overt act of
	relinquishment in perpetuity of all present and future rights to this
	software under copyright law.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
	IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
	OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
	OTHER DEALINGS IN THE SOFTWARE.

	For more information, please refer to <http://unlicense.org/>

	Richard Geldreich, Jr.
	12/30/2021
*/
/*
MIT License

Copyright (c) 2017  Joe Hood

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include<iostream>
#include<stdio.h>
#include<cmath>
#include<vector>

using namespace std;

namespace  AsciiPlotterUtils {

int max(vector<int> &data)
{
	int xmax = data[0];
    for (size_t i = 1; i < sizeof(data); i++)
	{
		if (data[i] > xmax)
		{
			xmax = data[i];
		}
	}
	return xmax;
}

int min(vector<int> &data)
{
	int xmin = data[0];
    for (size_t i = 1; i < data.size(); i++)
	{
		if (data[i] < xmin)
		{
			xmin = data[i];
		}
	}
	return xmin;
}

double max(vector<double> &data)
{
	double xmax = data[0];
    for (size_t i = 1; i < data.size(); i++)
	{
		if (data[i] > xmax)
		{
			xmax = data[i];
		}
	}
	return xmax;
}

double min(vector<double> &data)
{
	double xmin = data[0];
    for (size_t i = 1; i < data.size(); i++)
	{
		if (data[i] < xmin)
		{
			xmin = data[i];
		}
	}
	return xmin;
}

int max(int x1, int x2)
{
	if (x1 > x2)
	{
		return x1;
	}
	else
	{
		return x2;
	}
}

int min(int x1, int x2)
{
	if (x1 < x2)
	{
		return x1;
	}
	else
	{
		return x2;
	}
}

double max(double x1, double x2)
{
	if (x1 > x2)
	{
		return x1;
	}
	else
	{
		return x2;
	}
}

double min(double x1, double x2)
{
	if (x1 < x2)
	{
		return x1;
	}
	else
	{
		return x2;
	}
}

double map(double x, double in_min, double in_max, double out_min, double out_max)
{
	return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

vector<double> resample(vector<double> ydata, int newlength)
{
	int oldlength = ydata.size();
	vector<double> newdata(newlength);
	double factor = 1.0;
	double x, y, x1, y1, x2, y2;

	if (oldlength == newlength)
	{
		return ydata;
	}
	else
	{
		factor = (double)oldlength / (double)newlength;

		for (int newindex = 0; newindex < newlength; newindex++)
		{
			x = (double)newindex * factor;
			x1 = floor(x);
			x2 = x1 + 1.0;

			y1 = ydata[min(max(0, (int)x1), oldlength - 1)];
			y2 = ydata[min(max(0, (int)x2), oldlength - 1)];

			y = y1 + (y2 - y1) * (x - x1) / (x2 - x1);

			newdata[min(max(0, newindex), newlength - 1)] = y;
		}

		newdata[0] = ydata[0];
		newdata[newlength - 1] = ydata[oldlength - 1];

		return newdata;
	}
}

}

AsciiPlotter::AsciiPlotter()
{
	_title = "";
	_width = 100;
	_height = 50;
	_curves = 0;
}

AsciiPlotter::AsciiPlotter(string title)
{
	_title = title;
	_width = 100;
	_height = 50;
	_curves = 0;
}

AsciiPlotter::AsciiPlotter(string title, int width, int height)
{
	_title = title;
	_width = width;
	_height = height;
	_curves = 0;
}

void AsciiPlotter::addPlot( vector<double> ydata, string label = "", char marker = ' ')
{
	_markers[_curves] = marker;
	_labels[_curves] = label;
	_ydata[_curves++] = ydata;
}

void AsciiPlotter::show()
{
    double xmin = 0.0;
    double xmax = 1.0;

	double ymax = 1.0e-15;
	double ymin = 1.0e15;

	vector<double> resampled;
	vector< vector<char> > plane;

	int padding;
	string lmargin = "          ";

	plane.resize(_width);
	for (int i = 0; i < _width; i++)
	{
		plane[i].resize(_height);
		for (int j = 0; j < _height; j++)
		{
			plane[i][j] = ' ';
		}
	}

	for (int curve = 0; curve < _curves; curve++)
	{
        double mx = AsciiPlotterUtils::max(_ydata[curve]);
        double mn = AsciiPlotterUtils::min(_ydata[curve]);

		if (mx > ymax)
		{
			ymax = mx;
		}

		if (mn < ymin)
		{
			ymin = mn;
		}
	}

	for (int curve = 0; curve < _curves; curve++)
	{

        resampled = AsciiPlotterUtils::resample(_ydata[curve], _width);

		for (int row = 0; row < _width; row++)
		{
            int colindex = (int)AsciiPlotterUtils::map(resampled[row], ymin, ymax, 0.0, (double)_height);
			plane[row][min(max(0, colindex), _height - 1)] = _markers[curve];
		}
	}

	// title:

	cout << endl << endl;

    for (size_t i = 0; i < lmargin.length() + (_width - _title.length()) / 2 - 1; i++)
	{
		cout << " ";
	}

	cout << _title << endl << endl;
	
	// main plot plane:

	printf(" %8.2g +", ymax);

	for (int row = 0; row < _width; row++)
	{
		cout << "-";
	}
	cout << "+" << endl;

	for (int col = _height - 1; col >= 0; col--)
	{
		if (col == _height / 2)
		{
			padding = lmargin.length() - _ylabel.length();
			if (padding >= 0)
			{
				int totwidth = 0;
				for (int i = 0; i < padding - 1; i++)
				{
					cout << " ";
					totwidth++;
				}
				cout << _ylabel;
                for (size_t i = totwidth; i < lmargin.length() - _ylabel.length(); i++)
				{
					cout << " ";
				}
				cout << "|";

			}
			else
			{

			}
		}
		else
		{
			cout << lmargin << "|";
		}
		for (int row = 0; row < _width; row++)
		{
			cout << plane[row][col];
		}
		cout << "|" << endl;
	}

	printf(" %8.2g +", ymin);
	for (int row = 0; row < _width; row++)
	{
		cout << "-";
	}
	cout << "+" << endl;

	cout << lmargin;
	printf("%-8.2g", xmin);

	int buf = (_width - _xlabel.length()) / 2 - 6;
	for (int row = 0; row < buf; row++)
	{
		cout << " ";
	}

	cout << _xlabel;

	for (int row = 0; row < buf-1; row++)
	{
		cout << " ";
	}

	printf("%8.2g", xmax);

	cout << endl << endl;
	
	// legend:

	if (_legend)
	{
		cout << lmargin << "+";
		for (int row = 0; row < _width; row++)
		{
			cout << "-";
		}
		cout << "+" << endl;

		for (int curve = 0; curve < _curves; curve++)
		{
			cout << lmargin << "|   " << _markers[curve] << " " << _labels[curve];

            for (size_t i = 0; i < (_width - _labels[curve].length() - 5); i++)
			{
				cout << " ";
			}
			cout << "|" << endl;
		}

		cout << lmargin << "+";
		for (int row = 0; row < _width; row++)
		{
			cout << "-";
		}
		cout << "+" << endl;
	}
}

void AsciiPlotter::xlabel(string label)
{
	_xlabel = label;
}

void AsciiPlotter::ylabel(string label)
{
	_ylabel = label;
}

void AsciiPlotter::legend()
{
	_legend = true;
}

AsciiPlotter::~AsciiPlotter()
{
	// do nothing
}
