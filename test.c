#include <stdlib.h>
#include <stdio.h>
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

int peak_compare(const void * a, const void * b){ if(((peak_t *)b)->val > ((peak_t *)a)->val) return(1); else return(0); }

int main()
{

    const int Nsv = 34;
    const int OS = 4;
    const int Ndata = 65536; // length of complex baseband data record.
    const int Nsearch_fft = 16384;  // the size of the fft used for correlation must be at least twice as long as the data which is 8K.
    const float Fs = 4.0*1.023e6;
    const float thresh = 1600.0;
    FILE* fp;

    clock_t time0, time1;
    double exec_time;
    time0 = clock();


    ne10_fft_cpx_float32_t (*ca_float)[Nsearch_fft]=malloc(Nsv*Nsearch_fft*sizeof(ne10_fft_cpx_float32_t)); // ca_float[Nsv][Nsearch_fft]
    ne10_fft_cpx_float32_t (*CA_float)[Nsearch_fft]=malloc(Nsv*Nsearch_fft*sizeof(ne10_fft_cpx_float32_t)); // CA_float[Nsv][Nsearch_fft]
    uint8_t (*ca_seq)[1023*OS]=malloc(Nsv*1023*OS*sizeof(uint8_t ));                                        // ca_seq[Nsv][1023*OS]

    printf("computing C/A sequences.\n");
    for (int sv=1; sv<Nsv; sv++) cacode(sv, OS, ca_seq[sv]);

    //for (int i=0; i<10*OS; i++) printf("%+d,",ca_seq[1][i]); printf("\n");
    
    // convert to bpsk floats and zero pad
    for (int sv=1; sv<Nsv; sv++){
        for (int i=0; i<1023*OS; i++){     // convert to bpsk, 0 -> +1.0, 1 -> -1.0
            if (ca_seq[sv][i] == 0) {
                ca_float[sv][i].r = +1.0; 
                ca_float[sv][i].i =  0.0; 
            } else {
                ca_float[sv][i].r = -1.0;
                ca_float[sv][i].i =  0.0;
            }
        }
        for (int i=1023*OS; i<Nsearch_fft; i++){ // zero pad
            ca_float[sv][i].r = 0.0;
            ca_float[sv][i].i = 0.0;
        }
    }

    //for (int i=0; i<10*OS; i++) printf("%+2.1f,",ca_float[2][i].r); printf("\n");

    // Initialise Ne10, using hardware auto-detection to set library function pointers
    if (ne10_init() != NE10_OK) {
        fprintf(stderr, "Failed to initialise Ne10.\n");
        return 1;
    } else {
        printf("NE10_OK\n");
    }

    printf("computing fft twiddle factors\n");
    ne10_fft_cfg_float32_t search_fft_cfg;
    search_fft_cfg = ne10_fft_alloc_c2c_float32_c(Nsearch_fft);
    ne10_fft_cfg_float32_t long_fft_cfg;
    long_fft_cfg = ne10_fft_alloc_c2c_float32_c(Ndata);

    printf("computing C/A fft's\n");
    for (int sv=1; sv<Nsv; sv++) {
        ne10_fft_c2c_1d_float32_c(CA_float[sv], ca_float[sv], search_fft_cfg, 0);   // fft
        for (int i=0; i<Nsearch_fft; i++) CA_float[sv][i].i *= -1.0;                // complex conjugate.
    }

    //fp = fopen("CA_float.dat","w"); for (int i=0; i<Nsearch_fft; i++) fprintf(fp, "%f  %f\n", CA_float[1][i].r,CA_float[1][i].i); fclose(fp);

    // compute table of doppler frequencies
    const int Ndopp = 41;
    const float Dopp_start = -4000.0;
    const float Dopp_step  = 200.0;
    float dopp_index[Ndopp];
    printf("creating table of doppler frequencies\n");
    for (int i=0; i<Ndopp; i++) dopp_index[i] = Dopp_start + i*Dopp_step;

    //for (int i=0; i<Ndopp; i++) printf("%+f\n", dopp_index[i]);

    printf("computing the array of doppler mix vectors\n");
    ne10_fft_cpx_float32_t mix[Ndopp][Nsearch_fft/2];
    for (int k=0; k<Ndopp; k++) {
        //printf("computing mix vector for Doppler = %lf\n", dopp_index[k]);
        // compute the complex sinusiod at the negative of the doppler.
        float phase;
        for (int i=0; i<Nsearch_fft/2; i++){
            phase = -2.0*M_PI*(dopp_index[k]/Fs)*i;  // -2*pi*f*i
            mix[k][i].r = cos(phase);
            mix[k][i].i = sin(phase);
        }
    }

    //fp = fopen("mix.dat","w"); for (int i=0; i<Nsearch_fft/2; i++) fprintf(fp, "%f  %f\n", mix[25][i].r,mix[25][i].i); fclose(fp);

    time1 = clock();
    exec_time = (double)(time1-time0)/CLOCKS_PER_SEC;
    printf("initialization time = %lf\n", exec_time);


    // *******************
    //  Everything until this point could be done offline in advance.
    // *******************


    // read in the baseband data, 64K samples. 8k samples are used in search. 64k are used for refined doppler calculation.
    printf("reading the ADC data\n");
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
    peak_t peak_array[Nsv][Ndopp]; 
    printf("looping over the doppler bins.\n");
    for (int k=0; k<Ndopp; k++) {

        printf("\r%8.2lf", dopp_index[k] ); fflush(stdout);

        // compute the dot product to mix off the doppler, then zero pad.
        ne10_fft_cpx_float32_t s_mix[Nsearch_fft];
        for (int i=0; i<Nsearch_fft/2; i++)            {s_mix[i] = cpx_mul(s[i], mix[k][i]);}   // multiply
        for (int i=Nsearch_fft/2; i<Nsearch_fft; i++)  {s_mix[i].r = 0.0; s_mix[i].i = 0.0;}    // zero pad

        // convert mixed adc data to frequency domain.
        ne10_fft_cpx_float32_t S_mix[Nsearch_fft];
        ne10_fft_c2c_1d_float32_c(S_mix, s_mix, search_fft_cfg, 0);

        // loop over the SV computing the dot product of the frequency domain of the signal and C/A code. Then take the ifft() and find the peak.
        for (int sv=1; sv<Nsv; sv++) {

            // multiply the frequency domain representations of the signal and the ca code.
            ne10_fft_cpx_float32_t S_corr[Nsearch_fft];
            for (int i=0; i<Nsearch_fft; i++)   S_corr[i] = cpx_mul(S_mix[i], CA_float[sv][i]);  // complex multiply

            // take the inverse fft
            ne10_fft_cpx_float32_t s_corr[Nsearch_fft];
            ne10_fft_c2c_1d_float32_c(s_corr, S_corr, search_fft_cfg, 1); 

            // find the peak of the absolute value, save value and index.
            peak_array[sv][k] = sat_peak_find(s_corr, Nsearch_fft/2);

        }
    }

    fp = fopen("peak.dat","w"); 
    for (int sv=1; sv<Nsv; sv++) {
        for (int k=0; k<Ndopp; k++) {
            fprintf(fp, "%f ", peak_array[sv][k].val);
        } 
        fprintf(fp, "\n");
    }
    fclose(fp);

    // Now we have the correlation peak for each SV and each doppler bin.
    // For each SV wee need to find the doppler bin with the max correlation
    // and save a list of the SV that meet threshold.
    peak_t max_peak[Nsv];
    max_peak[0].val = 0.0; max_peak[0].loc = 0;
    int num_sat = 0;
    for (int sv=1; sv<Nsv; sv++) {
        max_peak[sv].val = 0.0; max_peak[sv].loc = 0; max_peak[sv].sv = sv;
        for (int k=0; k<Ndopp; k++) {
            if (peak_array[sv][k].val > max_peak[sv].val) {
                max_peak[sv].val = peak_array[sv][k].val;
                max_peak[sv].loc = peak_array[sv][k].loc;
            }
        }
        if (max_peak[sv].val > thresh) num_sat++;
    }

    time1 = clock();
    exec_time = (double)(time1-time0)/CLOCKS_PER_SEC;
    printf("\nsearch time = %lf\n", exec_time);

    qsort(max_peak, (size_t)Nsv, sizeof(peak_t), peak_compare);  // sort list of satellites by strength. the first num_sat entries are the SV found.

    printf("num_sat = %d\n", num_sat);
    for (int sv=0; sv<Nsv; sv++) if (max_peak[sv].val > thresh) printf("SV = %3d, max val = %5.2f, max loc = %5d\n", max_peak[sv].sv, max_peak[sv].val, max_peak[sv].loc); 


    // **********************
    // Now the list of visible satellites is known with their code delay offsets. We can use this to get a refined doppler estimate.
    // *********************

    // Refined doppler is found by first de-spreading the data with the code delay for each PRN in the list.
    // Then a 64K point fft can be computed to determine the doppler. 64K samples = 16ms. 1/16ms = 60Hz fft resolution so error should be in the range -30Hz to +30Hz.
    //
    // For each SV in the list, replicate the ca sequence to be longer than 64K plus the max shift value of 4092.
    // Then compute the dot product of the raw data (64K samples) with the portion of that long sequence corresponding to the code shift.
    // Then take the FFT and find the frequency bin with the highest absolute value. That bin corresponds to the exact doppler of the SV.
    printf("computing refined doppler estimate\n");
    const int Nrep = 2+Ndata/(1023*OS);
    int8_t ca_long_seq[Nrep*1023*OS];
    ne10_fft_cpx_float32_t s_despread[Ndata];
    for (int i=0; i<num_sat; i++){
        // fill in the first copy of the ca sequence
        for (int j=0; j<1023*OS; j++) if (ca_seq[max_peak[i].sv][j]==1) ca_long_seq[j] = -1; else ca_long_seq[j] = +1;
        // fill in the rest of the long sequence
        for (int k=1; k<Nrep; k++) for (int j=0; j<1023*OS; j++) ca_long_seq[k*1023*OS+j] = ca_long_seq[j];
        // de-spread
        for (int j=0; j<Ndata; j++) {
            s_despread[j].r = s_full[j].r*ca_long_seq[j+max_peak[i].loc];
            s_despread[j].i = s_full[j].i*ca_long_seq[j+max_peak[i].loc];
        }
        // calculate fft

        // find peak and convert to Hz
    }


    free(ca_float);
    free(ca_seq);
    free(CA_float);

}
