#ifndef ARSENAL_H
#define ARSENAL_H

using namespace std;

unsigned long int random_seed();
int sgn(double val);
void logspace(double * x, double a, double b, int n);
void linspace(double * x, double a, double b, int n);
void stratify_npts(double a, double b, int n1, int npts, double * x);
void scalepoints(double * x, double a, double b, double scale, int n);

//miscellaneous functions needed for interpolation routines
long binarySearch(double * A, int length, double value, bool skip_out_of_range);
void get_1D_derivatives(double * x, double * f, double * derivatives, int length, double);
void get_2D_derivatives(double * x, double * y, double ** f, double ** f1, double ** f2, double ** f12, int x_length, int y_length, double default_edge_fill);
void bcucof(double * y, double * y1, double * y2, double * y12, double d1, double d2, double ** c);
void bcuint(double * y, double * y1, double * y2, double * y12, double x1l, double x1u, double x2l, double x2u, double x1, double x2, double &ansy, double &ansy1, double &ansy2);
void polint(double xa[], double ya[], long n, double x, double *y, double *dy);
void polin2(double * x1a, double * x2a, double ** ya, long m, long n, double x1, double x2, double *y, double *dy);

//individual interpolation routines
double interpLinearDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpLinearNondirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiLinearDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiLinearNondirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpTriLinearDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpTriLinearNondirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpCubicDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
//double interpCubicNondirect(double * x, double * y, double xx, long size);
double interpCubicNonDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiCubicDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size);
double interpBiCubicDirectALT(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiCubicNonDirectALT(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpPolyDirect(double * x, double * y, double x0, long size);
double interpBiPolyDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size);

//miscellaneous interpolation routines
double interpNewtonDirect(double * x, double * y, double x0, long size);

//main interpolation routines
double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);
double interpolate2D(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);
double interpolate3D(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
			long x_size, long y_size, long z_size, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);

#endif
