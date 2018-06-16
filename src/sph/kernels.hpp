#include <cmath>

namespace kernels {

    struct Result2D {
       double w = 0;
       double dwdx = 0;
       double dwdy = 0;
    };

    struct Result1D {
       double w = 0;
       double dw = 0;
    };

    const float PI = 3.14159;

    double yang(double r, double h)  {
        double k = 2.0;
        if (fabs(r) < 2.0 * h) {
            return (4.0 * cos(PI/(2.0*h) * r) + cos(PI/h * r) + 3.0) * 1.0/(6.0*k*h);
        }
        return 0;
    }

    Result2D yang(double x, double y, double r, double h) {
        Result2D res;
        double q = fabs(r) / h;
        double k = 2.0;
        double sig = PI/((3*PI*PI-16)*(k*h)*(k*h));
        if (q <= 2.0) {
            res.w = (4 * cos(PI/2*q) + cos(PI * q) + 3) * sig;
            double derConst = (2 * sig * PI / (h * sqrt(x*x + y*y))) * (2*sin(PI)/2*q + sin(PI*q));
            res.dwdx = x * derConst;
            res.dwdy = y * derConst;
        }
        return res;
    }

    double wendlandC2(double r, double h) {
        double q = fabs(r) / h;
        if (q < 2.0) {
            return pow(1.0 - 0.5*q, 3.0) * (1.5*q + 1) * 5.0/(8.0*h);
        }
        return 0;
    }

    double wendlandC4(double r, double h) {
        double q = fabs(r) / h;
        if (q < 2.0) {
            return pow(1.0 - 0.5*q, 5.0) * (2.0*q*q + 2.5*q + 1) * 3.0/(4.0*h);
        }
        return 0;
    }

    double wendlandC6(double r, double h) {
        double q = fabs(r) / h;
        if (q < 2.0) {
            return pow(1.0 - 0.5*q, 7.0) * (21.0/8.0*q*q*q + 19.0/4.0*q*q + 3.5*q + 1) * 55.0/(64.0*h);
        }
        return 0;
    }

    Result1D splineCubic(double r, double h) {
        Result1D res;
        double q = fabs(r / h);
        double sig = 2/(3.0*h);
        if (q <= 1.0) {
            res.w = (1 - 1.5 * q*q + 0.75 * q*q*q) * sig;
            res.dw = (-3.0 * q * (1 - 0.5 * q) + 0.75 * q*q) * sig;
            if (r < 0) res.dw *= -1;
        } else if (q <= 2.0) {
            res.w = 0.25 * pow(2 - q, 3) * sig;
            res.dw = -0.75 * pow(2 - q, 2) * sig;
            if (r < 0) res.dw *= -1;
        }
        return res;
    }

    Result2D splineCubic(double x, double y, double r, double h) {
        Result2D res;
        double q = fabs(r) / h;
        double sig = 10/(7*PI*h*h);
        if (q <= 1.0) {
            res.w = (1 - 1.5 * q*q + 0.75 * q*q*q) * 2/(3.0*h);
            res.dwdx = (6*x*(0.5*q-1)+1.5*h*q*x) * (sig/h*h);
            res.dwdy = (6*y*(0.5*q-1)+1.5*h*q*y) * (sig/h*h);
        } else if (q <= 2.0) {
            res.w = (0.25 * (2 - q)*(2 - q)*(2 - q)) * 2/(3.0*h);
            res.dwdx = (0.5*q-1)*(0.5*q-1) * x * (-3*sig/h);
            res.dwdy = (0.5*q-1)*(0.5*q-1) * y * (-3*sig/h);
        }
        return res;
    }

    double splineQuintic(double r, double h) {
        double q = fabs(r) / h;
        if (q <= 1.0) {
            return (pow(3.0 - q, 5.0) - 6.0 * pow(2.0 - q, 5.0) + 15.0 * pow(1.0 - q, 5.0)) * 1/(120.0*h);
        } else if (q <= 2.0) {
            return (pow(3.0 - q, 5.0) - 6.0 * pow(2.0 - q, 5.0)) * 1/(120.0*h);
        } else if (q <= 3.0) {
            return pow(3.0 - q, 5.0) * 1/(120.0*h);
        }
        return 0;
    }


}
