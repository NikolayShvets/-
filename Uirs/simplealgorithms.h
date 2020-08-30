#ifndef SIMPLEALGORITHMS_H
#define SIMPLEALGORITHMS_H

#include <math.h>
#include "linearalgebra.h"

class simpleAlgorithms
{
public:
    static long double mZero(){
        long double temp{1}, result{1};
        while (result + 1.0 > 1) {
            result = temp;
            temp /= 2.0;
        }
        return result;
    }

    static long double getMax(const long double &a, const long double &b){
        if (a > b){
            return a;
        }else{
            return  b;
        }
    }

    static long double getMin(const long double &a, const long double &b){
        if (a < b){
            return a;
        }else{
            return  b;
        }
    }
};

#endif // SIMPLEALGORITHMS_H
