/*=============================================================================
   Copyright (c) 2021 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#pragma once
#ifndef SPLINE_H
#define SPLINE_H

#include<iostream>
#include<vector>
#include<unordered_map>
#include<cmath>
#include<functional>
#include<memory.h>
#include"utilities.h"

using namespace std;
namespace hypercurve {

template<typename T>
class cubic_spline
{
public:

    inline double *spl_memory_incr(size_t s)
    {
        double *ptr = _spl_mem_incr;
        _spl_mem_incr += s;
        return ptr;
    }

    void process(double *ptr, size_t size, memory_vector<control_point> &pts, double * spline_memory, size_t spline_mem_size)
    {

        interp.init(ptr, size, false);
        n = (pts.size() - 1);
        _spline_memory =  spline_memory;
        _spline_mem_size = (spline_mem_size);
        _spl_mem_incr = (spline_memory);
        x.init(spl_memory_incr(n + 1), n + 1, false);
        y.init(spl_memory_incr(n + 1), n + 1, false);
        a.init(spl_memory_incr(n), n, false);
        b.init(spl_memory_incr(n), n, false);
        c.init(spl_memory_incr(n), n, false);
        d.init(spl_memory_incr(n), n, false);
        h.init(spl_memory_incr(n), n, false);
        sig.init(spl_memory_incr(n + 1), n + 1, false);
        sig_temp.init(spl_memory_incr(n - 1), n - 1, false);
        tridiagonal.init( _spl_mem_incr, (n - 1) * n, false);

        interpolate_from_points(pts, size);

    }


    cubic_spline() {}
    cubic_spline(double *ptr, size_t size, memory_vector<control_point> &pts, double * spline_memory, size_t spline_mem_size)
        : interp(ptr, size, false)
        , n(pts.size() - 1)
        , _spline_memory(spline_memory)
        , _spline_mem_size(spline_mem_size)
        , _spl_mem_incr(spline_memory)
        , x(spl_memory_incr(n + 1), n + 1, false)
        , y(spl_memory_incr(n + 1), n + 1, false)
        , a(spl_memory_incr(n), n, false)
        , b(spl_memory_incr(n), n, false)
        , c(spl_memory_incr(n), n, false)
        , d(spl_memory_incr(n), n, false)
        , h(spl_memory_incr(n), n, false)
        , sig(spl_memory_incr(n + 1), n + 1, false)
        , sig_temp(spl_memory_incr(n - 1), n - 1, false)
        , tridiagonal( _spl_mem_incr, (n - 1) * n, false)
    {
        interpolate_from_points(pts, size);
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
                int iii = (n * i) + i;
                int kii = (n * k) + i;
                //double term = tridiagonal[k][i] / tridiagonal[i][i];
                double term = tridiagonal[kii] / tridiagonal[iii];
                for(j = 0; j < n; j++)
                {

                    int kk = n * k;
                    int jj = j;
                    int ii = n * i;
                    int fi = kk + jj;
                    int li = jj + ii;
                    tridiagonal[fi] = tridiagonal[fi] - term * tridiagonal[li];
                    /*
                    tridiagonal[k][j] = tridiagonal[k][j] -
                            term * tridiagonal[i][j];
                    */
                }
            }
        }

        // Back substitution
        for(i = m - 1; i >= 0; i--)
        {
            int in = (i * n) + (n - 1);
            sig_temp[i] = tridiagonal[in];
            //sig_temp[i] = tridiagonal[i][n - 1];
            for(j = i + 1; j < n - 1; j++)
            {
                int ij = (i * n) + j;
                //sig_temp[i] = sig_temp[i] - tridiagonal[i][j] * sig_temp[j];
                sig_temp[i] = sig_temp[i] - tridiagonal[ij] * sig_temp[j];
            }
            int ii = (i * n) + i;
            sig_temp[i] = sig_temp[i] / tridiagonal[ii];
            //sig_temp[i] = sig_temp[i] / tridiagonal[i][i];
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
            double ii = (i * n) + i;
            tridiagonal[ii] = 2 * (h[i] + h[i + 1]);
            //tridiagonal[i][i] = 2 * (h[i]+ h[i + 1]);
        }
        for(i = 0; i < n - 2; i++)
        {
            int ifirst = (n * i) + (i + 1);
            tridiagonal[ifirst] = h[i +  1];
            int ilast = (i+1) * n + i;
            tridiagonal[ilast] = h[i + 1];
            //tridiagonal[i][i + 1] = h[i + 1];
            //tridiagonal[i + 1][i] = h[i + 1];
        }
        for(i = 1; i < n; i++)
        {
            int ind = (i - 1) * n + (n - 1);
            tridiagonal[ind] = (y[i + 1] - y[i]) * 6.0 / h[i] -
                    (y[i] - y[i - 1]) * 6.0 / h[i - 1];
            /*
            tridiagonal[i - 1][n - 1] = (y[i + 1] - y[i]) * 6.0 / h[i] -
                    (y[i] - y[i - 1]) * 6.0 / h[i - 1];
            */
        }
    }

    T equation(T xval, T x_i, int i)
    {
        return (a[i] * pow((xval - x_i), 3)) + (b[i] * pow((xval - x_i), 2)) +
                (c[i] * (xval - x_i)) + d[i];
    }

    void calculate_interpolation()
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

    /*
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
    */

    template<typename N>
    void sn(N a)
    {
        std::cout << std::to_string(a) << std::endl;
    }

    void  reset(int n, memory_vector<point> &p)
    {
        /*
        sn(n);
        sn(x.size());
        sn(y.size());
        sn(p.size());
        */
        for(size_t  i = 0; i < p.size(); i++)
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

        ::memset(tridiagonal.data(), 0, n * (n-1) * sizeof(T));

        ::memset(interp.data(), 0, sizeof(T) * n_precision);
        /*
        for(int i = 0; i < tridiagonal.size(); i++)
        {
            ::memset(tridiagonal[i].data(), 0, n * sizeof(T) );
        }
        */

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
        calculate_interpolation();
    }

    void interpolate_from_points(memory_vector<point> &p,
                                            int precision)
    {
        int n = p.size() - 1;
        if(n_points != p.size())
        {
            //this->resize(n);
            n_points = p.size();
        }
        if(precision != n_precision)
        {
            n_precision = precision;
        }

        reset(n, p);

        this->interpolate(n);
    }

protected:
    int n_points = 0, n_precision = 0;
    int n;

    memory_vector<double> interp;
    double *_spline_memory;
    double *_spl_mem_incr;
    size_t _spline_mem_size;

    memory_vector<T> x, y, a, b, c, d, sig, sig_temp, h;
    memory_vector<T> tridiagonal;
};



}


#endif // SPLINE_H
