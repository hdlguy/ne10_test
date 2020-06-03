
struct sat_peak_struct {ne10_float32_t val; int loc; int sv;};
typedef struct sat_peak_struct peak_t;

double sat_abs(ne10_fft_cpx_float32_t num);
struct sat_peak_struct sat_peak_find(ne10_fft_cpx_float32_t *vec, int N);


