#ifndef FUNCTIONS_H
#define FUNCTIONS_H


#include <complex>
#include <math.h>
#include "params.h"
#include <boost/math/special_functions/bessel.hpp>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_specfunc.h>

using namespace std;
using namespace boost::math;

complex<double> I = 1i;

inline const auto H (const auto& x) ->
    typename std::remove_reference<decltype(x)>::type {
    return x > 0 ? 1 : 0;
}

class IntegratorQAGI{
public:
    IntegratorQAGI(){
        ws = gsl_integration_workspace_alloc(1000);
    }
    ~IntegratorQAGI(){
        gsl_integration_workspace_free(ws);
    }
    double integrate(double f(double,void*), void* p, double sing){
        func.function = f;
        func.params = p;
        gsl_integration_qag(&func, 0, sing - eps, 1e-6, 1e-6, 1000, GSL_INTEG_GAUSS15, ws,&result1, &error);
        //cout<<error<<endl;
        gsl_integration_qagiu(&func, sing + eps, 1e-6, 1e-6, 1000, ws, &result2, &error);
        //cout<<error<<endl;
        return result1 + result2;
    }

private:
    gsl_integration_workspace *ws;
    gsl_function func;
    double result1;
    double result2;
    double error;

};

struct sum_params{
    double psi;
    double theta;
    double r;
    double mu;
    int m;
};

double sum_elem (int m, double s, double psi, double theta, double r, double mu){
    return -(4 * m * pow(s,2) * sin(m * (psi - theta)) * omega/
            (pow(m, 2)*pow(omega, 2)-pow(s, 2)*pow(U, 2))*
               (
                (mu>r)?
                gsl_sf_bessel_In_scaled(m,r*abs(s))*gsl_sf_bessel_Kn_scaled(m,mu*abs(s))*exp(-abs(s)*(mu-r)):
                gsl_sf_bessel_Kn_scaled(m,r*abs(s))*gsl_sf_bessel_In_scaled(m,mu*abs(s))*exp(-abs(s)*(r-mu))
                ));
}

double undint_s(double s, void *p){
    sum_params *par = (sum_params* )p;
    return sum_elem(par->m, s ,par->psi, par->theta, par->r, par->mu);
}

IntegratorQAGI qagi;

double integr_s(int m, double psi, double theta, double r, double mu){
    sum_params par;
    par.psi = psi;
    par.theta = theta;
    par.r = r;
    par.mu = mu;
    par.m = m;
    return qagi.integrate(&undint_s, &par, m*omega/U);
}



double sum_integr (double psi, double theta, double r, double mu){
    double sprev = 0;
    double sm = 0;
    int m = 1;
    const int m_max = 30;
    do {
        sprev = sm;
        sm = sm + integr_s (m, psi, theta, r, mu);
        ++m;
    } while (m < m_max && abs ((sm - sprev)/sprev)> err);
    return sm;
}

#endif // FUNCTIONS_H
