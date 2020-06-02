
const int Nsv = 34;
const int OS = 4;
const int Ndata = 65536; // length of complex baseband data record.
const int Nsearch_fft = 16384;  // the size of the fft used for correlation must be at least twice as long as the data which is 8K.

struct sat_peak_struct {double val; int loc;};

double sat_abs(ne10_fft_cpx_float32_t num);
struct sat_peak_struct sat_peak_find(ne10_fft_cpx_float32_t *vec, int N);


