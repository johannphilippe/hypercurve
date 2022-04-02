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
#include"utilities.h"

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

    void calculate_interpolation(int n)
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
                                            int precision, point size)
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
