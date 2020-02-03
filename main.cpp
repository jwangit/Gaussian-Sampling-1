/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 * 
 */

#include <cstdlib>
#include <math.h>
#include <NTL/RR.h>
#include "Samplers.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}//end-get_timestamp()

int main(void) {
    
    timestamp_t ts_start, ts_end;
    RR sigma = to_RR(3);
    RR c = to_RR(0);
    int precision = 64;
    int tailcut = 11;

    
    Samplers *samplers = new Samplers(sigma, c, tailcut, precision);
    samplers->BuildProbabilityMatrix();
    timestamp_t avgUsual = 0.0;
    
    int nIterations = 10;//, tailcut;
 //   long precision;

    // (Roy, Vercauteren and Verbauwhede, 2013): precision and sigma values for statistical distance less than 2^{-90}
//    precision = 107;
//   tailcut = 13;
            
/*
 * Usual Gaussian sampler
 */
                        
	int sample=0;
	nIterations = 1;    
        for(int it = 0; it < nIterations; it++) {
           ts_start = get_timestamp();
	     
           sample = samplers->KnuthYao();
           ts_end = get_timestamp();

           avgUsual += (ts_end - ts_start);

         //  cout << "[!] KuthYao  sampler running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
         //  cout << "[>] Sample from the lattice: " << sample << endl;

            }//end-for


            if(nIterations > 1)
                cout << "[!] KnuthYao Gaussian sampler average running time: " << (float)(((float)(nIterations)*1000000000.0)/avgUsual) << " samples/s.\n" << endl;
           

    delete(samplers);
    
    return 0;
    
}//end-main() 
