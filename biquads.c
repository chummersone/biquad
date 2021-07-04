#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "biquads.h"

#define BQN (3)

// index type
typedef unsigned int index_t;
// coeff type
typedef double coeff_t;

// the biquad structure
struct biquad_s {
    coeff_t   A[BQN];
    coeff_t   B[BQN];
    sample_t  X[BQN];
    sample_t  Y[BQN];
    index_t   index;
    double    fs;
    bq_type_e type;
    double    fc;
    double    Q;
    double    peakGain;
};

// private functions
static double complex biquad_evalPoly(double complex z, const coeff_t coeffs[BQN]);
static double complex biquad_tf(const biquad* b, double freq);
static void biquad_calculate(biquad* b);

// return buffer index
static inline index_t buff_ix(index_t n, int k) {
    return (n + k + BQN) % BQN;
}

// create a biquad
biquad* biquad_create(
    bq_type_e type,
    int fs,
    double fc,
    double Q,
    double peakGain
) {
    // allocate biquad and return pointer
    biquad* b = malloc(sizeof(biquad));
    if (b != NULL) {
    
        // put data in structure
        b->fs = fs;
        b->index = 0;
        b->type = type;
        b->fc = fc;
        b->Q = Q;
        b->peakGain = peakGain;
        
        // clear buffers
        biquad_clear(b);
        
        // calculate filter
        biquad_calculate(b);
    }
    return b;
}

// destroy the biquad
void biquad_destroy(biquad* b) {
    free(b);
}

// clone the biquad
biquad* biquad_clone(const biquad* original) {
    biquad* b = malloc(sizeof(biquad));
    if (b != NULL) {
        *b = *original;
    }
    return b;
}

// calculate magnitude
double biquad_magnitude(const biquad* b, double freq) {
    return cabs(biquad_tf(b, freq));
}

// calculate wrapped phase
double biquad_wrphase(const biquad* b, double freq) {
    double phase = carg(biquad_tf(b, freq));
    return (phase < 0 ? phase + (2.0*M_PI) : phase);
}

// process sample
void biquad_process(biquad* b, sample_t* buffer, unsigned int num_samples) {

    coeff_t *A  = b->A;
    coeff_t *B  = b->B;
    sample_t *X = b->X;
    sample_t *Y = b->Y;
    index_t nn = b->index;

    for (unsigned int i = 0; i < num_samples; i++) {
        index_t nm1 = buff_ix(nn, -1);
        index_t nm2 = buff_ix(nn, -2);

        // put input on to buffer
        X[nn] = buffer[i];

        // process input
        Y[nn] =
            B[0] * X[nn] +
            B[1] * X[nm1] +
            B[2] * X[nm2] -
            A[1] * Y[nm1] -
            A[2] * Y[nm2];

        // write output
        buffer[i] = Y[nn];

        nn = buff_ix(nn, 1);
    }

    b->index = nn;
}

// set filter parameter - fs
void biquad_setFs(biquad* b, double fs) {
    b->fs = fs;
    biquad_calculate(b);
}

// set filter parameter - type
void biquad_setType(biquad* b, bq_type_e type) {
    b->type = type;
    biquad_calculate(b);
}

// set filter parameter - fc
void biquad_setFc(biquad* b, double fc) {
    b->fc = fc;
    biquad_calculate(b);
}

// set filter parameter - Q
void biquad_setQ(biquad* b, double Q) {
    b->Q = Q;
    biquad_calculate(b);
}

// set filter parameter - peakGain
void biquad_setPeakGain(biquad* b, double peakGain) {
    b->peakGain = peakGain;
    biquad_calculate(b);
}

// get filter parameter - fs
double biquad_getFs(const biquad* b) {
    return b->fs;
}

// get filter parameter - type
bq_type_e biquad_getType(const biquad* b) {
    return b->type;
}

// get filter parameter - fc
double biquad_getFc(const biquad* b) {
    return b->fc;
}

// get filter parameter - Q
double biquad_getQ(const biquad* b) {
    return b->Q;
}

// get filter parameter - peakGain
double biquad_getPeakGain(const biquad* b) {
    return b->peakGain;
}

// clear buffers
void biquad_clear(biquad* b) {
    // reset buffers
    for (int n = 0; n < BQN; n++) {
        b->X[n] = 0.0;
        b->Y[n] = 0.0;
    }
}

// evaluate complex polynomial
double complex biquad_evalPoly(double complex z, const coeff_t coeffs[BQN]) {
    double complex y = 0;
    for (int i = 0; i < BQN; i++) {
        y += (coeffs[i]) * cpow(z, i);
    }
    return y;
}

// calculate transfer function
double complex biquad_tf(const biquad* b, double freq) {

    /* Calculate z = e^{\omega T} for frequency */
    double complex omegaT = 2 * M_PI * freq / (b->fs);
    double complex z = cexp(I*omegaT);
    
    /* Calculate the value of the numerator in the transfer function */
    double complex num = biquad_evalPoly(z,b->B);
    
    /* Calculate the value of the denominator in the transfer function */
    double complex den = biquad_evalPoly(z,b->A);

    return num / den;
}

// calculate filter coefficients
void biquad_calculate(biquad* b) {

    // shorten names
    coeff_t *A      = b->A;
    coeff_t *B      = b->B;
    double Q        = b->Q;
    double fs       = b->fs;
    double fc       = b->fc;
    double peakGain = b->peakGain;

    coeff_t AA       = pow(10.0, peakGain / 40.0);
    coeff_t w0       = 2.0 * M_PI * fc / fs;
    coeff_t alpha    = sin(w0) / (2.0 * Q);
    
    coeff_t cos_w0   = cos(w0);
    coeff_t sqrt_AA  = sqrt(AA);
    
    // source : http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
    
    switch (b->type) {
        case BQ_LOWPASS:
            B[0] = (1.0 - cos_w0) / 2.0;
            B[1] = 1.0 - cos_w0;
            B[2] = (1.0 - cos_w0) / 2.0;
            A[0] = 1 + alpha;
            A[1] = -2.0 * cos_w0;
            A[2] = 1.0 - alpha;
            break;
        case BQ_HIGHPASS:
            B[0] = (1.0 + cos_w0) / 2.0;
            B[1] = -(1.0 + cos_w0);
            B[2] = (1.0 + cos_w0) / 2.0;
            A[0] = 1.0 + alpha;
            A[1] = -2.0 * cos_w0;
            A[2] = 1.0 - alpha;
            break;
        case BQ_BANDPASS: // (constant 0 dB peak gain)
            B[0] = alpha;
            B[1] = 0.0;
            B[2] = -alpha;
            A[0] = 1.0 + alpha;
            A[1] = -2.0 * cos_w0;
            A[2] = 1.0 - alpha;
            break;
        case BQ_NOTCH:
            B[0] = 1.0;
            B[1] = -2.0 * cos_w0;
            B[2] = 1.0;
            A[0] = 1.0 + alpha;
            A[1] = -2.0 * cos_w0;
            A[2] = 1.0 - alpha;
            break;
        case BQ_PEAK:
            B[0] = 1.0 + alpha*AA;
            B[1] = -2.0 * cos_w0;
            B[2] = 1.0 - alpha*AA;
            A[0] = 1.0 + alpha/AA;
            A[1] = -2.0 * cos_w0;
            A[2] = 1.0 - alpha/AA;
            break;
        case BQ_LOWSHELF:
            B[0] = AA * ( (AA + 1.0) - (AA - 1.0) * cos_w0 + 2.0 * sqrt_AA * alpha);
            B[1] = 2.0 * AA * ((AA - 1) - (AA + 1.0) * cos_w0);
            B[2] = AA * ((AA + 1.0) - (AA - 1.0) * cos_w0 - 2.0 * sqrt_AA * alpha);
            A[0] = (AA + 1.0) + (AA - 1.0) * cos_w0 + 2.0 * sqrt_AA * alpha;
            A[1] = -2.0 * ((AA - 1.0) + (AA + 1.0) * cos_w0);
            A[2] = (AA + 1.0) + (AA - 1.0) * cos_w0 - 2.0 * sqrt_AA * alpha;
            break;
        case BQ_HIGHSHELF:
            B[0] = AA * ((AA + 1.0) + (AA - 1.0) * cos_w0 + 2.0 * sqrt_AA * alpha);
            B[1] = -2.0 * AA * ((AA - 1.0) + (AA + 1.0) * cos_w0);
            B[2] = AA * ((AA + 1.0) + (AA - 1.0) * cos_w0 - 2.0 * sqrt_AA * alpha);
            A[0] = (AA + 1.0) - (AA - 1.0) * cos_w0 + 2.0 * sqrt_AA * alpha;
            A[1] = 2.0 * ((AA - 1.0) - (AA + 1.0) * cos_w0);
            A[2] = (AA + 1.0) - (AA - 1.0) * cos_w0 - 2.0 * sqrt_AA * alpha;
            break;
        case BQ_NONE:
        case BQ_NUM_FILTERS:
        default:
            B[0] = 1.0;
            B[1] = 0.0;
            B[2] = 0.0;
            A[0] = 1.0;
            A[1] = 0.0;
            A[2] = 0.0;
    }
    
    // normalize
    coeff_t norm = A[0];
    for (index_t i = 0; i < BQN; i++) {
        A[i] /= norm;
        B[i] /= norm;
    }
}

// return the filter name
const char* biquad_getFilterName(const biquad* b) {
    switch (b->type) {
        case BQ_LOWPASS:
            return "low-pass";
        case BQ_HIGHPASS:
            return "high-pass";
        case BQ_BANDPASS: // (constant 0 dB peak gain)
            return "band-pass";
        case BQ_NOTCH:
            return "notch";
        case BQ_PEAK:
            return "peak";
        case BQ_LOWSHELF:
            return "low-shelf";
        case BQ_HIGHSHELF:
            return "high-shelf";
        case BQ_NONE:
        case BQ_NUM_FILTERS:
        default:
            return "none";
    }
}
