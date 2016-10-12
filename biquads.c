#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "biquads.h"

#define BQN 3

// return buffer index
#define BUFFIX(ix,n) ((ix+n+BQN)%BQN)

// the biquad structure
struct biquad_s {
    double A[BQN];
    double B[BQN];
    bqFloat X[BQN];
    bqFloat Y[BQN];
    unsigned int index;
    double fs;
    bqFilterType type;
    double fc;
    double Q;
    double peakGain;
};

// private functions
static double complex biquad_evalPoly(double complex z, double coeffs[]);
static double complex biquad_tf(biquad *b, double freq);
static void biquad_calculate(biquad *b);

// create a biquad
biquad *biquad_create(
    bqFilterType type,
    int fs,
    double fc,
    double Q,
    double peakGain
) {
    // allocate biquad and return pointer
    biquad *b = malloc(sizeof(biquad));
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
void biquad_destroy(biquad *b) {
    free(b);
}

// calculate magnitude
double biquad_magnitude(biquad *b, double freq) {
    return cabs(biquad_tf(b, freq));
}

// process sample
bqFloat biquad_process(biquad *b, bqFloat input) {

    double *A = b->A;
    double *B = b->B;
    bqFloat *X = b->X;
    bqFloat *Y = b->Y;
    unsigned int ix = b->index;
    
    // put input on to buffer
    X[BUFFIX(ix,0)] = input;
    
    // process input
    Y[BUFFIX(ix,0)] =
        B[0]*X[BUFFIX(ix,  0)] +
        B[1]*X[BUFFIX(ix, -1)] +
        B[2]*X[BUFFIX(ix, -2)] -
        A[1]*Y[BUFFIX(ix, -1)] -
        A[2]*Y[BUFFIX(ix, -2)];
    
    // write output
    bqFloat output = Y[BUFFIX(ix,0)];
    
    // step through buffer
    ix = BUFFIX(ix,1);
    
    return output;
}

// set filter parameter - fs
void biquad_setFs(biquad *b, double fs) {
    b->fs = fs;
    biquad_calculate(b);
}

// set filter parameter - type
void biquad_setType(biquad *b, bqFilterType type) {
    b->type = type;
    biquad_calculate(b);
}

// set filter parameter - fc
void biquad_setFc(biquad *b, double fc) {
    b->fc = fc;
    biquad_calculate(b);
}

// set filter parameter - Q
void biquad_setQ(biquad *b, double Q) {
    b->Q = Q;
    biquad_calculate(b);
}

// set filter parameter - peakGain
void biquad_setPeakGain(biquad *b, double peakGain) {
    b->peakGain = peakGain;
    biquad_calculate(b);
}

// clear buffers
void biquad_clear(biquad *b) {
    // reset buffers
    for (int n = 0; n < BQN; n++) {
        b->X[n] = 0.0;
        b->Y[n] = 0.0;
    }
}

// evaluate complex polynomial
double complex biquad_evalPoly(double complex z, double coeffs[]) {
    double complex y = 0;
    for (int i = 0; i < BQN; i++) {
        y += (coeffs[i]) * cpow(z,i);
    }
    return y;
}

// calculate transfer function
double complex biquad_tf(biquad *b, double freq) {

    /* Calculate z = e^{\omega T} for frequency */
    double complex omegaT = 2 * M_PI * freq / (b->fs);
    double complex z = cexp(I*omegaT);
    
    /* Calculate the value of the numerator in the transfer function */
    double complex num = biquad_evalPoly(z,b->B);
    
    /* Calculate the value of the denominator in the transfer function */
    double complex den = biquad_evalPoly(z,b->A);
    
    return num/den;
}

// calculate filter coefficients
void biquad_calculate(biquad *b) {

    // shorten names
    double *A = b->A;
    double *B = b->B;
    double Q = b->Q;
    double fs = b->fs;
    double fc = b->fc;
    double peakGain = b->peakGain;

    double AA  = pow(10.0,peakGain/40.0);
    double w0 = 2.0*M_PI*fc/fs;
    double alpha = sin(w0)/(2.0*Q);
    
    // source : http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
    
    switch (b->type) {
        case BQ_LOWPASS:
            B[0] = (1.0 - cos(w0))/2.0;
            B[1] = 1.0 - cos(w0);
            B[2] = (1.0 - cos(w0))/2.0;
            A[0] = 1 + alpha;
            A[1] = -2.0*cos(w0);
            A[2] = 1.0 - alpha;
            break;
        case BQ_HIGHPASS:
            B[0] = (1.0 + cos(w0))/2.0;
            B[1] = -(1.0 + cos(w0));
            B[2] = (1.0 + cos(w0))/2.0;
            A[0] = 1.0 + alpha;
            A[1] = -2.0*cos(w0);
            A[2] = 1.0 - alpha;
            break;
        case BQ_BANDPASS: // (constant 0 dB peak gain)
            B[0] = alpha;
            B[1] = 0.0;
            B[2] = -alpha;
            A[0] = 1.0 + alpha;
            A[1] = -2.0*cos(w0);
            A[2] = 1.0 - alpha;
            break;
        case BQ_NOTCH:
            B[0] = 1.0;
            B[1] = -2.0*cos(w0);
            B[2] = 1.0;
            A[0] = 1.0 + alpha;
            A[1] = -2.0*cos(w0);
            A[2] = 1.0 - alpha;
            break;
        case BQ_PEAK:
            B[0] = 1.0 + alpha*AA;
            B[1] = -2.0*cos(w0);
            B[2] = 1.0 - alpha*AA;
            A[0] = 1.0 + alpha/AA;
            A[1] = -2.0*cos(w0);
            A[2] = 1.0 - alpha/AA;
            break;
        case BQ_LOWSHELF:
            B[0] = AA*( (AA+1.0) - (AA-1.0)*cos(w0) + 2.0*sqrt(AA)*alpha );
            B[1] = 2.0*AA*( (AA-1) - (AA+1.0)*cos(w0) );
            B[2] = AA*( (AA+1.0) - (AA-1.0)*cos(w0) - 2.0*sqrt(AA)*alpha );
            A[0] = (AA+1.0) + (AA-1.0)*cos(w0) + 2.0*sqrt(AA)*alpha;
            A[1] = -2.0*( (AA-1.0) + (AA+1.0)*cos(w0) );
            A[2] = (AA+1.0) + (AA-1.0)*cos(w0) - 2.0*sqrt(AA)*alpha;
            break;
        case BQ_HIGHSHELF:
            B[0] = AA*( (AA+1.0) + (AA-1.0)*cos(w0) + 2.0*sqrt(AA)*alpha );
            B[1] = -2.0*AA*( (AA-1.0) + (AA+1.0)*cos(w0) );
            B[2] = AA*( (AA+1.0) + (AA-1.0)*cos(w0) - 2.0*sqrt(AA)*alpha );
            A[0] = (AA+1.0) - (AA-1.0)*cos(w0) + 2.0*sqrt(AA)*alpha;
            A[1] = 2.0*( (AA-1.0) - (AA+1.0)*cos(w0) );
            A[2] = (AA+1.0) - (AA-1.0)*cos(w0) - 2.0*sqrt(AA)*alpha;
            break;
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
    double norm = A[0];
    for (int i = 0; i < BQN; i++) {
        A[i] = A[i]/norm;
        B[i] = B[i]/norm;
    }
}
