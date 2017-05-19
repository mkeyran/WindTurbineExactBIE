#ifndef PARAMS_H
#define PARAMS_H

#include <array>
#include <utility>

const int n_sect_y = 50;
const int n_sect_z = 12;

const int N2 = n_sect_y;
const int N3 = n_sect_z;

const double eps = 1e-6;

struct geom_sect {
    double b_begin;
    double b_end;
    double aust;
};

using geom_t = std::array<geom_sect, n_sect_y>;
geom_t geom;

geom_t&& make_geom(){
    geom_t tarr;
    for (int i=0;i<n_sect_y; ++i){
        tarr[i].b_begin = -0.138; //Координата задней кромки сечения
        tarr[i].b_end = 0.138; // Координата передней кромки сечения
        tarr[i].aust = 4. * M_PI / 180; // Угол установки
    }
    return std::move(tarr);
}

struct length {
    double begin;
    double end;
} l {0.08, 0.8};

const double omega = 1400. * 2. * M_PI / 60;
const double U = 8;

const double tol = 1e-20;
const double err = 1e-6;

#endif // PARAMS_H
