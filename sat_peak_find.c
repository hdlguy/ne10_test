#include <math.h>
#include "NE10.h"
#include "sat_search.h"

struct sat_peak_struct sat_peak_find (ne10_fft_cpx_float32_t x[], int N)
{

    struct sat_peak_struct peak;
    peak.loc = 0; peak.val = 0.0;
    ne10_float32_t peak_val = 0;
    int peak_loc = 0;
    for (int32_t shift=0; shift<N; shift++) {
        ne10_float32_t val = x[shift].r*x[shift].r + x[shift].i*x[shift].i;
        if (val>peak_val) {
            peak_val = val;
            peak_loc = shift;
        }
    }

    peak.val = sqrt(peak_val); // take square root once at end.
    peak.loc = peak_loc;

    return(peak);

}


