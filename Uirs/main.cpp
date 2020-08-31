#include <iostream>
#include "custommodel.h"
#include "integrator.h"
#include <fstream>

using namespace std;

int main()
{
    satellite *s = new satellite();
    s->setSampIncrement(0.1);
    s->setT_st(0.0);
    s->setT_fin(30.0);
    dormandPrinceIntgrator *dp_integrator = new dormandPrinceIntgrator();
    eulerIntegrator *e = new eulerIntegrator();
    dp_integrator->run(s);
    //e->run(s);
    print(s->getResult());
    return 0;
}
