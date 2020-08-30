#ifndef IINTQRATOR_H
#define IINTQRATOR_H

#include "mathmodel.h"
#include "simplealgorithms.h"
#include "fstream"
#include <vector>

class IIntegrator
{
protected:
    long double eps;
public:
    IIntegrator() : eps( 1e-16 ) {}
    inline void setPrecision( long double eps ) { this->eps = eps; }
    inline long double getPrecision() const { return eps; }
    virtual void run(IMathModel *model) = 0;
};

class eulerIntegrator : public IIntegrator
{
public:
    eulerIntegrator();
    void run(IMathModel *model) override;
};

class dormandPrinceIntgrator : public IIntegrator
{
private:
    vector c, b1, b2;
    matrix a;
    std::vector<vector> K;
    long double zero;
public:
    dormandPrinceIntgrator();
    void run(IMathModel *model);
};

#endif // IINTQRATOR_H
