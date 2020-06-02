#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "NE10.h"
#include <math.h>
#include <time.h> 
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

    const int Nsv = 34;
    const int OS = 4;
    const int Ndata = 65536; // length of complex baseband data record.
    const int Nsearch_fft = 16384;  // the size of the fft used for correlation must be at least twice as long as the data which is 8K.

    clock_t time0, time1;
    double exec_time;
    time0 = clock();


    ne10_fft_cpx_float32_t (*ca_float)[Nsearch_fft]=malloc(Nsv*Nsearch_fft*sizeof(ne10_fft_cpx_float32_t)); // ca_float[Nsv][Nsearch_fft]
    ne10_fft_cpx_float32_t (*CA_float)[Nsearch_fft]=malloc(Nsv*Nsearch_fft*sizeof(ne10_fft_cpx_float32_t)); // CA_float[Nsv][Nsearch_fft]
    uint8_t (*ca_seq)[1023*OS]=malloc(Nsv*1023*OS*sizeof(uint8_t ));                                        // ca_seq[Nsv][1023*OS]

    printf("computing C/A sequences.\n");
    for (int sv=1; sv<Nsv; sv++) cacode(sv, OS, ca_seq[sv]);

    // convert to bpsk floats
    for (int sv=1; sv<Nsv; sv++){
        for (int i=0; i<1023*OS; i++){
            if (ca_seq[sv][i] == 0) {
                ca_float[sv][i].r = -1.0; 
                ca_float[sv][i].i =  0.0; 
            } else {
                ca_float[sv][i].r = +1.0;
                ca_float[sv][i].i =  0.0;
            }
        }
        for (int i=1023*OS; i<Nsearch_fft; i++){
            ca_float[sv][i].r = 0.0;
            ca_float[sv][i].i = 0.0;
        }
    }

    // Initialise Ne10, using hardware auto-detection to set library function pointers
    if (ne10_init() != NE10_OK) {
        fprintf(stderr, "Failed to initialise Ne10.\n");
        return 1;
    } else {
        printf("NE10_OK\n");
    }

    printf("configuring search fft\n");
    ne10_fft_cfg_float32_t search_fft_cfg;
    search_fft_cfg = ne10_fft_alloc_c2c_float32_c(Nsearch_fft);

    printf("computing C/A fft's\n");
    for (int sv=1; sv<Nsv; sv++) {
        ne10_fft_c2c_1d_float32_c(CA_float[sv], ca_float[sv], search_fft_cfg, 0);   // fft
        for (int i=0; i<Nsearch_fft; i++) CA_float[sv][i].i *= -1.0;                // complex conjugate.
    }

    // compute table of doppler frequencies
    const int Ndopp = 41;
    const float Dopp_start = -4000.0;
    const float Dopp_step  = 200.0;
    float dopp_index[Ndopp];
    printf("creating table of doppler frequencies\n");
    for (int i=0; i<Ndopp; i++) dopp_index[i] = Dopp_start + i*Dopp_step;

    printf("computing the array of doppler mix vectors\n");
    ne10_fft_cpx_float32_t mix[Ndopp][Nsearch_fft/2];
    for (int k=0; k<Ndopp; k++) {
        //printf("computing mix vector for Doppler = %lf\n", dopp_index[k]);
        // compute the complex sinusiod at the negative of the doppler.
        float phase;
        for (int i=0; i<Nsearch_fft/2; i++){
            phase = -2.0*M_PI*dopp_index[k]*i;  // -2*pi*f*i
            mix[k][i].r = cos(phase);
            mix[k][i].i = sin(phase);
        }
    }

    time1 = clock();
    exec_time = (double)(time1-time0)/CLOCKS_PER_SEC;
    printf("initialization time = %lf\n", exec_time);


    // *******************
    //  Everything until this point could be done offline in advance.
    // *******************


    // read in the baseband data, 64K samples. 8k samples are used in search. 64k are used for refined doppler calculation.
    printf("reading the ADC data\n");
    FILE* fp;
    fp = fopen("dataout.txt","r");
    ne10_fft_cpx_float32_t s_full[Ndata]; // holds the full length complex baseband data for doppler refinement.
    for (int i=0; i<Ndata; i++) {
        int rsamp, isamp;
        fscanf(fp, "%d, %d\n", &rsamp, &isamp);
        s_full[i].r = (ne10_float32_t) (rsamp);
        s_full[i].i = (ne10_float32_t) (isamp);
    }
    fclose(fp);

    time0 = clock();

    ne10_fft_cpx_float32_t s[Nsearch_fft/2]; // holds the 8K length complex baseband data for search.
    for (int i=0; i<Nsearch_fft/2; i++) {
        s[i] = s_full[i];
    }

    // loop over the doppler table computing correlations
    struct sat_peak_struct peak_array[Nsv][Ndopp]; 
    printf("looping over the doppler bins.\n");
    for (int k=0; k<Ndopp; k++) {

        printf("\r%8.2lf", dopp_index[k] ); fflush(stdout);

        // compute the dot product to wipe of the doppler and zero pad.
        ne10_fft_cpx_float32_t s_mix[Nsearch_fft];
        for (int i=0; i<Nsearch_fft/2; i++){
            s_mix[i] = cpx_mul(s[i], mix[k][i]);  // multiply
        }
        for (int i=Nsearch_fft/2; i<Nsearch_fft; i++){
            s_mix[i].r = 0.0; s_mix[i].i = 0.0;  // zero pad
        }

        // convert mixed adc data to frequency domain.
        ne10_fft_cpx_float32_t S_mix[Nsearch_fft];
        ne10_fft_c2c_1d_float32_c(S_mix, s_mix, search_fft_cfg, 0);

        // loop over the SV computing the dot product of the frequency domain of the signal and C/A code. Then take the ifft().
        for (int sv=1; sv<Nsv; sv++) {

            // multiply the frequency domain representations of the signal and the ca code.
            ne10_fft_cpx_float32_t S_corr[Nsearch_fft];
            for (int i=0; i<Nsearch_fft; i++){
                S_corr[i] = cpx_mul(S_mix[i], CA_float[k][i]);  // complex multiply
            }

            // take the inverse fft
            ne10_fft_cpx_float32_t s_corr[Nsearch_fft];
            ne10_fft_c2c_1d_float32_c(s_corr, S_corr, search_fft_cfg, 1); 

            // find the peak of the absolute value, save value and index.
            peak_array[sv][k] = sat_peak_find(s_corr, Nsearch_fft/2);

        }
    }

    // Now we have the correlation peak for each SV and each doppler bin.
    // For each SV wee need to find the doppler bin with the max correlation
    // and save a list of the SV that meet threshold.
    struct sat_peak_struct max_peak[Nsv];
    for (int sv=1; sv<Nsv; sv++) {
        max_peak[sv].val = 0.0;
        max_peak[sv].loc = 0;
        for (int k=0; k<Ndopp; k++) {
            if (peak_array[sv][k].val > max_peak[sv].val) {
                max_peak[sv].val = peak_array[sv][k].val;
                max_peak[sv].loc = peak_array[sv][k].loc;
            }
        }
    }

    time1 = clock();
    exec_time = (double)(time1-time0)/CLOCKS_PER_SEC;
    printf("\nsearch time = %lf\n", exec_time);

    for (int sv=1; sv<Nsv; sv++) {
        printf("SV=%3d, max val = %5.2f, max loc = %5d\n", sv, max_peak[sv].val, max_peak[sv].loc); 
    }

    free(ca_float);
    free(ca_seq);
    free(CA_float);

}
