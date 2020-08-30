#include "custommodel.h"

satellite::satellite():IMathModel(){
    X0.resize(6);
    X0[0] = R + 0.0;    //X0
    X0[1] = R + 700.0;  //Y0
    X0[2] = R + 0.0;    //Z0
    X0[3] = 80.0;       //Vx0
    X0[4] = 20.0;       //Vy0
    X0[5] = 0.0;        //Vz0

    perturbation = new matrix(3, 3);
    //матрица имеет следующий вид: нулевая строка - встречный ветер, первая - попутный, вторая - боковой
    perturbation->operator()(0, 0) = -10.0; perturbation->operator()(0, 1) = 0.0; perturbation->operator()(0, 2) = 0.0;
    perturbation->operator()(1, 0) = 10.0; perturbation->operator()(1, 1) = 0.0; perturbation->operator()(1, 2) = 0.0;
    perturbation->operator()(2, 0) = 10.0; perturbation->operator()(2, 1) = 0.0; perturbation->operator()(2, 2) = -10.0;
}

void satellite::getRP(const vector &X,long double t, vector &Y) const{
    Y.resize(X.size());
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    //Рассчитываем расстояние от центра Земли до объекта
    long double ro = sqrt(pow(X[0],2.) + pow(X[1], 2.) + pow(X[2], 2.));
    //три слагаемых - соответствующие компоненты встречного, попутного и бокового ветров
    Y[3] = -nu*m*X[0]/pow(ro,3.) + perturbation->operator()(0, 0) + perturbation->operator()(1, 0) + perturbation->operator()(2, 0);
    Y[4] = -nu*m*X[1]/pow(ro,3.) + perturbation->operator()(0, 1) + perturbation->operator()(1, 1) + perturbation->operator()(2, 1);
    Y[5] = -nu*m*X[2]/pow(ro,3.) + perturbation->operator()(0, 2) + perturbation->operator()(1, 2) + perturbation->operator()(2, 2);
}

