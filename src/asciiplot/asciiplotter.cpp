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
#include"asciiplotter.h"

using namespace std;

int max(vector<int> data)
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

int min(vector<int> data)
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

double max(vector<double> data)
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

double min(vector<double> data)
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
		double mx = max(_ydata[curve]);
		double mn = min(_ydata[curve]);

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

		resampled = resample(_ydata[curve], _width);

		for (int row = 0; row < _width; row++)
		{
			int colindex = (int)map(resampled[row], ymin, ymax, 0.0, (double)_height);
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
