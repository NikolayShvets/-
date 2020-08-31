#include "integrator.h"
#define max(a, b) ( ( (a) > (b) ) ? (a) : (b) )
#define min(a, b) ( ( (a) < (b) ) ? (a) : (b) )



eulerIntegrator::eulerIntegrator():IIntegrator(){

}

void eulerIntegrator::run(IMathModel *model)
{
    std::ofstream file;
    file.open("results.txt");
    model->clearResult();
    model->prepareResult();
    long double
            t = model->getT_st(),
            t_out = model->getT_st(),
            t1 = model->getT_fin(),
            h = model->getSampIncrement();
    vector
            X = model->getInitialConditions(),
            Y(X.size());
    while(t < t1){
        model->getRP(X, t, Y);
        Y = X + Y*h;
        X = Y;
        if (t_out < t + h) model->addResult(X, t);
        t += h;
    }
}

const long double dormandPrinceIntgrator::c[7] = { 0, 1./5, 3./10, 4./5, 8./9, 1., 1. };
const long double dormandPrinceIntgrator::a[7][6] = {
    { 0. },
    { 1./5 },
    { 3./40, 9./40 },
    { 44./45, -56./15, 32./9 },
    { 19372./6561, -25360./2187, 64448./6561, -212./729 },
    { 9017./3168, -355./33, 46732./5247, 49./176, -5103./18656 },
    { 35./384, 0., 500./1113, 125./192, -2187./6784, 11./84 }
};
const long double dormandPrinceIntgrator::b1[7] = { 35./384, 0., 500./1113, 125./192, -2187./6784, 11./84, 0 };
const long double dormandPrinceIntgrator::b2[7] = { 5179./57600, 0., 7571./16695, 393./640, -92097./339200, 187./2100, 1./40 };


dormandPrinceIntgrator::dormandPrinceIntgrator(): IIntegrator(){
    zero = simpleAlgorithms::mZero();
}

void dormandPrinceIntgrator::run(IMathModel * model)
{
    /*std::ofstream file;
    file.open("results.txt");
    long double
            t = model->getT_st(),
            t_out = t,
            t1 = model->getT_fin(),
            h,
            h_new = model->getSampIncrement(),
            e = 0;

    vector
            X = model->getInitialConditions(),
            X1( X.size() ),
            X2( X.size() ),
            Xout ( X.size() ),
            Y( X.size() );
    model->clearResult();
    model->prepareResult();
    for(auto &elem : K){
        elem.resize(X.size());
    }

    while ( t < t1 )
    {
        h = h_new;
        for (int j = 0; j < (int)K.size(); ++j) {
            for (int k = 0; k < X.size(); ++k) {
                Y[k] = X[k];
                for (int s = 0; s < j; ++s) {
                    Y[k] += K[s][k] * a(j,s) * h;
                }
            }
            model->getRP(Y, t, K[j]);
        }
        e = 0;
        for (int i = 0; i < X.size(); ++i){
            X1[i] = X2[i] = X[i];
            for (int j = 0; j < b1.size(); ++j) {
                X1[i] += K[j][i] * b1[j] * h;
                X2[i] += K[j][i] * b2[j] * h;
            }
            e += powl(h * (X1[i] - X2[i]) /
                      simpleAlgorithms::getMax(simpleAlgorithms::getMax(abs(X[i]), fabsl(X1[i])), simpleAlgorithms::getMax(1e-5L, 2 * zero / eps)), 2.0);
        }
        e = sqrtl( e / X.size() );
        h_new = h / simpleAlgorithms::getMax( 0.1, simpleAlgorithms::getMin( 5., pow(e / eps, 0.2)/0.9 ) );

        if ( e > eps )
            continue;
        while ( (t_out < t + h) && (t_out <= t1) )
        {
            long double theta = (t_out - t)/h, b[6];

            b[0] = theta  * ( 1 + theta    *(-1337./480.  + theta*(1039./360.    + theta*(-1163./1152.))));
            b[1] = 0;
            b[2] = 100.   * powl(theta, 2) * (1054./9275. + theta*(-4682./27825. + theta*(379./5565.)))/3.;
            b[3] = -5.    * powl(theta, 2) * (27./40.     + theta*(-9./5.        + theta*(83./96.)))/2.;
            b[4] = 18225. * powl(theta, 2) * (-3./250.    + theta*(22./375.      + theta*(-37./600.)))/848.;
            b[5] = -22.   * powl(theta, 2) * (-3./10.     + theta*(29./30.       + theta*(-17./24.)))/7.;

            for ( int k = X.size()-1; k >= 0; k-- )
            {
                long double sum  = 0;
                for ( int j = 5; j >= 0; j-- )
                    sum += b[j] * K[j][k];
                Xout[k] = X[k] + h * sum;
            }

            model->addResult( Xout, t_out );
            for(int i = 0; i < Xout.size(); i++){
                file<<Xout[i]<<"|";
            }
            file<<t_out<<std::endl;
            t_out += model->getSampIncrement();
        }
        X = X1;
        t += h;
    }*/
    std::ofstream file;
    file.open("results.txt");
    long double
                t = model->getT_st(),
                t_out = t,
                t1 = model->getT_fin(),
                h,
                h_new = model->getSampIncrement(),
                e = 0;

    vector
            X = model->getInitialConditions(),
            X1( X.size() ),
            X2( X.size() ),
            Xout ( X.size() ),
            Y( X.size() );
    model->prepareResult();
    for ( int j = 7; j > 0; --j, K[j].resize( X.size() ) );
    int N = 0;
    while ( t < t1 )
    {
        h = h_new;
        for ( int j = 0; j < 7; j++ )
        {
            for ( int k = X.size()-1; k >= 0; k-- )
            {
                Y[k] = X[k];
                for ( int s = 0; s < j; s++ )
                {
                    Y[k] += K[s][k] * a[j][s] * h;
                }
            }
           model->getRP( Y, t + c[j] * h, K[j] );
        }
        e = 0;
        for ( int k = X.size()-1; k >= 0; k-- )
        {
            X1[k] = X2[k] = X[k];
            for ( int j = 0; j < 7; j++ )
            {
                X1[k] += K[j][k] * b1[j] * h;
                X2[k] += K[j][k] * b2[j] * h;
            }
            e += pow( h * (X1[k] - X2[k]) / max( max( fabsl(X[k]), fabsl(X1[k]) ), max((long double)1e-5, 2*u/eps) ) , 2 );
        }
        e = sqrtl( e / X.size() );
        h_new = h / max( 0.1, min( 5., pow(e / eps, 0.2)/0.9 ) );

        if ( e > eps )
            continue;

        while ( (t_out < t + h) && (t_out <= t1) )
        {
            long double l_ldTheta = (t_out - t)/h,
                        b[6];

            b[0] = l_ldTheta * ( 1 + l_ldTheta*(-1337./480. + l_ldTheta*(1039./360. + l_ldTheta*(-1163./1152.))));
            b[1] = 0;
            b[2] = 100. * pow(l_ldTheta, 2) * (1054./9275. + l_ldTheta*(-4682./27825. + l_ldTheta*(379./5565.)))/3.;
            b[3] = -5. * pow(l_ldTheta, 2) * (27./40. + l_ldTheta*(-9./5. + l_ldTheta*(83./96.)))/2.;
            b[4] = 18225. * pow(l_ldTheta, 2) * (-3./250. + l_ldTheta*(22./375. + l_ldTheta*(-37./600.)))/848.;
            b[5] = -22. * pow(l_ldTheta, 2) * (-3./10. + l_ldTheta*(29./30. + l_ldTheta*(-17./24.)))/7.;

            for ( int k = X.size()-1; k >= 0; k-- )
            {
                long double l_ldSum  = 0;
                for ( int j = 5; j >= 0; j-- )
                    l_ldSum += b[j] * K[j][k];
                Xout[k] = X[k] + h * l_ldSum;
            }

            model->addResult( Xout, t_out );
            for(int i = 0; i < Xout.size(); i++){
                file<<Xout[i]<<"|";
            }
            file<<t_out<<std::endl;
            t_out += model->getSampIncrement();
        }
        X = X1;
        t += h;
        N++;
    }
}
