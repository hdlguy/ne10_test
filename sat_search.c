

double sat_abs(ne10_fft_cpx_float32_t num) {
    return( sqrt(num.i*num.i + num.r*num.r) );
}


// a function to find peak absolute value in a complex arrray.
sat_peak_struct sat_peak_find(ne10_fft_cpx_float32_t *vec, int N) {
    sat_peak_struct peak;
    peak.val = 0.0;
    peak.loc = 0;
    for (int i=0; i<N; i++) {
        double temp_val = sat_abs(vec[i]);
        if (temp_val > peak.val) {
            peak.val = temp_val;
            peak.loc = i;
        }
    }
    return(peak);
}

