/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <math.h>
#include "HIBE.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

int main(void) {
        
    int action = 2;
    
    timestamp_t ts_start, ts_end;
    int h, k;
    double lambda;
    int m1, m2, q, r;

    h = 10;
    k = 2; // n = 2^k is the degree of polynomials in R and R_0
    q = 11; //q must be prime and congruent to 3 mod 8
    m1 = 5;
    m2 = 30; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
    r = ceil((double)(1 + (log(q)/log(3))));
    lambda = ceil(1 + (log(q)/log(2)));

    if(q % 8 != 3) {
        cout << "q must be congruent to 3 mod 8.\n";
        return -1;
    }

    if(/*m1 < sigma ||*/ m1 < r) {
        cout << "m1 must be greater or equal to sigma and r.\n";
        return -1;        
    }

    if(m1 < lambda) {
        cout << "m1 must be greater or equal to lambda.\n";
        return -1;
    }

    if(m2 < lambda) {
        cout << "m2 must be greater or equal to lambda.\n";
        return -1;
    }

    if(m2 < lambda*m1) {
        cout << "m2 must be greater or equal to lambda times m1.\n";
        return -1;
    }

    HIBE *hibe = new HIBE(q, m1, m2, k);    
    
    RR sigmaRR = to_RR(3.195); // Standard deviation
    RR c = to_RR(0); // Center of the distribution            
    int tailcut = 13;
    int precision = 107;
    int nRectangles = 63; // Parameter of Ziggurat algorithm
    int omega = precision; // Parameter of Ziggurat algorithm
    
    switch(action) {
        case 1: {
            /* Parameter set from (Roy et al., 2013). That depends on the cryptographic system requirements */
            
            if(sigmaRR*to_RR(tailcut) > power2_RR(sizeof(int)*8+1)-1) {
                cout << "Error! This distribution can not be simulated. Aborting..." << endl;
                return -1;
            }//end-if                
            
            Vec<int> ZigguratPoly, KnuthPoly;
            int nSamples = 8194; // #coefficients in the polynomial
                        
            ts_start = get_timestamp();            
            ZigguratPoly = hibe->GetSampler()->PolyGeneratorZiggurat(nSamples, nRectangles, sigmaRR, omega, precision, tailcut); // Coefficients, rectangles, sigma, omega and precision
            ts_end = get_timestamp();            
                        
            cout << "[!] Ziggurat running time for " << nSamples << " samples: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << ZigguratPoly << endl;
            
            ts_start = get_timestamp();                        
            KnuthPoly = hibe->GetSampler()->PolyGeneratorKnuthYao(nSamples, precision, tailcut, sigmaRR, c); // Coefficients, precision, tailcut, and sigma
            ts_end = get_timestamp();            
            
            cout << "[!] Knuth-Yao running time for " << nSamples << " samples: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << KnuthPoly << endl;
            
            break;
        }
        case 2: {
            
            mat_RR BTilde;
            Vec<ZZX> A; // Block isometric basis A
            ZZX c, sample;
            RR normBTilde;
            ZZ zero;  
            
            c.SetMaxLength(hibe->GetN());
            zero = to_ZZ(0);
            
            for(int i = 0; i < hibe->GetN(); i++) // Center of distribution D_{\Lambda, \sigma, c}
                c[i] = zero;            
            
            /* Basis generation phase - (Mochetti, and Dahab, 2014) and (Stehlé et al., 2009) */
            ts_start = get_timestamp();
            hibe->Setup(h);
            ts_end = get_timestamp();            

            cout << "\n[!] Setup running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            
            cout << "\n/* Vector \\hat{a} generated by Setup algorithm */" << endl;
            cout << hibe->GetA() << endl;
            
            /* Rot_f(\hat{a}) - Expansion of vector \hat{a} into a matrix \*/
            hibe->GetSampler()->Rot(A, hibe->GetA(), hibe->GetM(), hibe->GetN());            
            
            cout << "\n/* Expansion of vector \\hat{a} into A */" << endl;            
            cout << A << endl;
            
            /* Orthogonalization of block isometric basis A - (Lyubashevsky, and Prest, 2015), Algorithm 5 */
            ts_start = get_timestamp();
            normBTilde = hibe->GetSampler()->BlockGSO(BTilde, A, hibe->GetM(), hibe->GetN(), precision);
            ts_end = get_timestamp();            

            cout << "\n[!] BlockGSO running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            
            // Uncomment this block if you wanna a preview of the Gram-Schmidt reduced basis
            /*cout << "\n[!] Norm of Gram-Schmidt reduced basis: " << normBTilde << endl;            
            for(int i = 0; i < hibe->GetM()*hibe->GetN(); i++) {
                cout << BTilde[i];
                if(i != 0 && i % hibe->GetM())
                    cout << endl;
            }//end-for*/
            
            sigmaRR = normBTilde*(log(hibe->GetN())/log(2)) + 1;
            
            /* Sampling from the discrete Gaussian distribution D_{\Lambda, \sigma, c} - (Lyubashevsky, and Prest, 2015), Algorithm 8 */
            ts_start = get_timestamp();    
            sample = hibe->GetSampler()->GaussianSamplerFromLattice(A, BTilde, sigmaRR, precision, tailcut, c, hibe->GetN());            
            ts_end = get_timestamp();    
            
            cout << "\n[!] GaussianSampler running time: " << (float)((ts_end - ts_start)/1000000.0) << " ms." << endl;
            cout << "\nSample from the lattice: " << sample << endl;
                        
            break;
        }
        default:
            break;            
    }//end-switch
    
    delete(hibe);
    return 0;
    
}