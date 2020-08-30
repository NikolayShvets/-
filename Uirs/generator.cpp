#include "generator.h"
using namespace std;
generator::generator()
{
    t_corr = 2.0L*M_PI/w;
}

long double generator::white_noise_generator(long double m, long double d)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<long double> distribution(m, d);
    nu_0 = distribution(generator);//*pow(D, 0.5);
    return nu_0;
}
generator::~generator()
{

}
