#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h> 
#include "NE10.h"
#include "sat_search.h"

void cacode(int sv, int OS, uint8_t* ca_up);

ne10_fft_cpx_float32_t cpx_mul(ne10_fft_cpx_float32_t a, ne10_fft_cpx_float32_t b){
    ne10_fft_cpx_float32_t temp;
    temp.r = a.r*b.r-a.i*b.i;
    temp.i = a.r*b.i+a.i*b.r;
    return(temp);
}


int main()
{

    // compute the array of ca_sequences at 4Msps
    ne10_fft_cpx_float32_t ca_float[Nsv][Nsearch_fft], CA_float[Nsv][Nsearch_fft]; // probably too big for the stack, 4092*34*8 = 1MB.
    uint8_t ca_seq[Nsv][1023*OS];
    printf("computing C/A sequences.\n");
    for (int sv=1; sv<Nsv; sv++){
        cacode(sv, OS, ca_seq[sv]);
    }

    for (int sv=1; sv<Nsv; sv++){
        for (int i=0; i<1023*OS; i++){
            if (ca_seq[sv][i] == 0) {
                ca_float[Nsv][i].r = -1.0; 
                ca_float[Nsv][i].i =  0.0; 
            } else {
                ca_float[Nsv][i].r = +1.0;
                ca_float[Nsv][i].i =  0.0;
            }
        }
        for (int i=1023*OS; i<Nsearch_fft; i++){
            ca_float[Nsv][i].r = 0.0;
            ca_float[Nsv][i].i = 0.0;
        }
    }

/*
    if (ne10_init() != NE10_OK) { fprintf(stderr, "Failed to initialise Ne10.\n"); return 1; }

    ne10_fft_cfg_float32_t search_fft_cfg;
    search_fft_cfg = ne10_fft_alloc_c2c_float32_c(Nsearch_fft);
    
    for (int sv=1; sv<Nsv; sv++){
        ne10_fft_c2c_1d_float32_c(CA_float[sv], ca_float[sv], search_fft_cfg, 0);
    }
*/

    // now we have the correlation templates converted to frequency domain.  All this could be done offline in advance.

/*
    // read in the baseband data, 64K samples. 8k samples are used in search. 64k are used for refined doppler calculation.
    // read the adc data.
    FILE* fp;
    fp = fopen("dataout.txt","r");
    ne10_fft_cpx_float32_t s_full[Ndata]; // holds the full length complex baseband data.
    for (int i=0; i<Ndata; i++) {
        int rsamp, isamp;
        fscanf(fp, "%d, %d\n", &rsamp, &isamp);
        s_full[i].r = (ne10_float32_t) (rsamp);
        s_full[i].i = (ne10_float32_t) (isamp);
    }
    fclose(fp);

    ne10_fft_cpx_float32_t s[Nsearch_fft/2]; // holds the 8K length complex baseband data.
    for (int i=0; i<Nsearch_fft/2; i++) {
        s[i] = s_full[i];
    }
*/
    
/*
    // compute table of doppler frequencies
    const int Ndopp = 41;
    const double Dopp_start = -4000.0;
    const double Dopp_step  = 200.0;
    double dopp_index[Ndopp];
    for (int i=0; i<Ndopp; i++) dopp_index[i] = Dopp_start + i*Dopp_step;
    //for (int i=0; i<Ndopp; i++) cout << dopp_index[i] << endl;
    
    // loop over the doppler table computing correlations
    printf("looping over the doppler bins.\n");
    for (int k=0; k<Ndopp; k++) {
        // compute the complex sinusiod at the negative of the doppler.
        double phase;
        ne10_fft_cpx_float32_t mix[Nsearch_fft/2];
        for (int i=0; i<Nsearch_fft/2; i++){
            phase = -2.0*M_PI*dopp_index[k]*i;  // -2*pi*f*i
            mix[i].r = cos(phase);
            mix[i].i = sin(phase);
        }

        // compute the dot product to wipe of the doppler and zero pad.
        ne10_fft_cpx_float32_t s_mix[Nsearch_fft];
        for (int i=0; i<Nsearch_fft/2; i++){
            s_mix[i] = cpx_mul(s[i], mix[i]);  // multiply
        }
        for (int i=Nsearch_fft/2; i<Nsearch_fft; i++){
            s_mix[i].r = 0.0; s_mix[i].i = 0.0;  // zero pad
        }

        // convert to frequency domain.

        // loop over the SV computing the dot product of the frequency domain of the signal and C/A code. Then take the ifft().
        for (int SV=1; SV<Nsv; SV++) {
            // multiply the frequency domain representations of the signal and the ca code.

            // take the inverse fft

            // find the peak of the absolute value, save value and index.

        }
    }
*/

}

