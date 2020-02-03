/* 
 * File:   Samplers.h
 * Author: jnortiz
 *
 * Created on April 24, 2015, 3:51 PM
 */

#include <NTL/mat_ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h>

#ifndef SAMPLERS_H
#define	SAMPLERS_H

using namespace std;
using namespace NTL;

class Samplers {
public:
       
    Samplers(RR sigma, RR c, int tailcut, int precision);
    virtual ~Samplers();
    
    /* Sampling from a discrete Gaussian distribution over the integers */
    int KnuthYao();
    void BuildProbabilityMatrix();
    RR Probability(RR x);

private:
        
    /* Knuth-Yao attributes */
    Vec< Vec<int> > P;
    Vec<int> begin;
    RR sigma;
    RR c;
    int tailcut;
    int precision;
    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    
};

#endif	/* SAMPLERS_H */

