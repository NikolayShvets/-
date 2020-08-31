#ifndef CUSTOMMODEL_H
#define CUSTOMMODEL_H

#include "mathmodel.h"
#include "simplealgorithms.h"

class satellite : public IMathModel
{
private:
    const long double nu{398600.436e9}; //гравитационный параметр Земли [km^3 * s^-1]
    long double m {2.0};                //масса тела [kg]
    long double r{0.4};                 //радиус тела [m]
    long double R{6371000.0};             //радус Земли [m]
    matrix *perturbation;               //матрица возмущений (ветер) 3х3 [m/s]
public:
    satellite();
    void getRP(const vector &X,long double t, vector &Y) const override;
};

#endif // CUSTOMMODEL_H
