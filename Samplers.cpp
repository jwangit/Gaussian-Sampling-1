/* 
 * File:   Samplers.cpp
 * Author: jnortiz
 * 
 * Created on April 24, 2015, 3:51 PM
 */

#include "Samplers.h"
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/RR.h>
#include <NTL/matrix.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <NTL/ZZX.h>

using namespace NTL;
using namespace std;

Samplers::Samplers(RR sigma, RR c, int tailcut, int precision){
     
this->sigma = sigma;
this->c = c;
this->tailcut = tailcut;
this->precision = precision;    
RR::SetPrecision(precision);
}//end-Samplers()

Samplers::~Samplers() {
    this->P.kill();
    this->begin.kill();
//    this->X.kill();
//    this->f.kill();
}

RR Samplers::Probability(RR x) {
    
    RR S = sigma*sqrt(2*ComputePi_RR());
    RR overS = 1/S;
    
    if(x == to_RR(0))
        return overS;
    
    return overS*exp(-(power((x-c)/sigma, 2))/2.0);
    
}//end-Probability()

// If the bit is zero, the output becomes "a" 
int Select(int a, int b, unsigned bit) {
    
    unsigned mask;
    int output;
    
    mask = -bit;
    output = mask & (a ^ b);
    output = output ^ a;
    
    return output;
    
}//end-Select()

int Samplers::KnuthYao() {

    int bound, center, col, d, invalidSample, pNumRows, pNumCols, S, signal;
    unsigned enable, hit;
    unsigned long r;
    
    bound = tailcut*to_int(sigma);
    center = to_int(c);
    d = 0;
    hit = 0;
    signal = 1 - 2*RandomBits_long(1);
    invalidSample = bound+1;
    pNumRows = this->P.length();
    pNumCols = this->P[0].length();    
    
    Vec<int> randomBits;
    randomBits.SetLength(pNumRows);
    
    int i, index, j, length;
    length = sizeof(unsigned long)*8; 
    
    index = 0;
    for(i = 0; i < (pNumRows/length+1); i++) {
        r = RandomWord();
        for(j = 0; j < length, index < pNumRows; j++, r >>= 1)
            randomBits[index++] = (r & 1);
    }//end-for
    
    S = 0;    
    for(int row = 0; row < pNumRows; row++) {
        
        d = 2*d + randomBits[row]; // Distance calculus
        
        for(col = this->begin[row]; col < pNumCols; col++) {
            
            d = d - this->P[row][col];
            
            enable = (unsigned)(d + 1); // "enable" turns 0 iff d = -1
            enable = 1 ^ ((enable | -enable) >> 31) & 1; // "enable" turns 1 iff "enable" was 0
             
            /* When enable&!hit becomes 1, "col" is added to "S";
             * e.g. enable = 1 and hit = 0 */
            S += Select(invalidSample, col, (enable & !hit));
            hit += (enable & !hit);
                            
        }//end-for
        
    }//end-for
    
    /* Note: the "col" value is in [0, bound]. So, the invalid sample must be 
     * greater than bound. */
    S %= invalidSample;
    S = S - bound + center;
    S *= signal;
    
    return S;
    
}//end-Knuth-Yao()

/* This method build the probability matrix for samples in the range 
 * [-tailcut*\floor(sigma), +tailcut*\floor(sigma)] */
void Samplers::BuildProbabilityMatrix() {
    
    RR::SetPrecision(to_long(precision));

    Vec< Vec<int> > auxP;
    Vec<int> auxBegin;
    
    // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
    int i, j, bound, pNumCols, pNumRows, x;
    vec_RR probOfX;
    RR pow;
    
    bound = tailcut*to_int(sigma);
    
    probOfX.SetLength(bound+1);
       
    auxP.SetLength(precision);
    for(i = 0; i < auxP.length(); i++)
        auxP[i].SetLength(bound+1);

    for(x = bound; x > 0; x--)
        probOfX[bound-x] = Probability(to_RR(x) + c);
    div(probOfX[bound], Probability(to_RR(0) + c), to_RR(2));
    
    i = -1;
    for(j = 0; j < precision; j++) {
        pow = power2_RR(i--); // 2^{i}
        for(x = bound; x >= 0; x--) {
            auxP[j][bound-x] = 0;                
            if(probOfX[bound-x] >= pow) {
                auxP[j][bound-x] = 1;
                probOfX[bound-x] -= pow;
            }//end-if
        }//end-for
    }//end-while
    
    this->P = auxP;
    
    // Uncomment this line if you want to preview the probability matrix P
    this->PrintMatrix("Probability matrix", this->P);
    
    pNumCols = this->P[0].length();
    pNumRows = this->P.length();
    
    auxBegin.SetLength(pNumRows);
    
    // Computing in which position the non-zero values in P start and end 
    for(i = 0; i < pNumRows; i++) {
        
        auxBegin[i] = pNumCols-1;
        
        for(j = 0; j < pNumCols; j++)
            if(this->P[i][j] == 1) {
                auxBegin[i] = j;
                break;
            }//end-if
        
    }//end-for
    
    this->begin = auxBegin;
                
}//end-BuildProbabilityMatrix()



/*vec_RR Samplers::GaussianSamplerFromLattice(const mat_ZZ& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, const vec_RR center) {

    cout << "\n[*] Gaussian Sampler status: ";
    
    RR::SetPrecision(precision);    

    double sizeof_RR = pow(2.0, sizeof(RR));
    RR d, innerp, innerp1, sigma_i, Z;
    vec_RR C, sample;    
    int cols, i, j, rows;
    
    cols = BTilde.NumCols();
    rows = BTilde.NumRows();
    
    C.SetLength(cols);
    sample.SetLength(cols);
                
    C = center;
    
    for(i = rows-1; i >= 0; i--) {
        
        NTL::InnerProduct(innerp, BTilde[i], BTilde[i]);
        NTL::InnerProduct(innerp1, C, BTilde[i]);        
        div(d, innerp1, innerp);        
        div(sigma_i, sigma, sqrt(innerp));
        
        if(floor(sigma_i) == 0)
            Z = d;
        else {
            if(sigma_i > sizeof_RR)
                sigma_i = to_RR(2)*sigma;
            this->BuildProbabilityMatrix(precision, tailcut, sigma_i, d);                     
            Z = to_RR(this->KnuthYao(tailcut, sigma_i, d));            
        }//end-else
        
        for(j = 0; j < B.NumCols(); j++)
            C[j] = C[j] - to_RR(B[i][j])*Z;
        
    }//end-for
    
    sub(sample, center, C);    
    cout << "Pass!" << endl;
    
    return sample;
    
}//end-GaussianSamplerFromLattice()
*/
void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()
