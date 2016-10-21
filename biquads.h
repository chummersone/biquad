#ifndef BIQUADS_H_
#define BIQUADS_H_

#ifdef __cplusplus
extern "C" {
#endif

// biquad type
typedef struct biquad_s biquad;

// biquad sample type (floating point)
typedef float bqFloat;

// filter type enumeration
typedef enum {
    BQ_NONE,
    BQ_LOWPASS,
    BQ_HIGHPASS,
    BQ_BANDPASS,
    BQ_NOTCH,
    BQ_PEAK,
    BQ_LOWSHELF,
    BQ_HIGHSHELF,
    BQ_NUM_FILTERS
} bqFilterType;

// create a biquad
extern biquad* biquad_create(
    bqFilterType type,
    int fs,
    double fc,
    double Q,
    double peakGain
);

// destroy the biquad
extern void biquad_destroy(biquad* b);

// clone the biquad
extern biquad* biquad_clone(const biquad* b);

// calculate magnitude
extern double biquad_magnitude(const biquad* b, double freq);

// process sample
extern bqFloat biquad_process(biquad* b, bqFloat input);

// set filter parameters
extern void biquad_setFs(biquad* b, double fs);
extern void biquad_setType(biquad* b, bqFilterType type);
extern void biquad_setFc(biquad* b, double fc);
extern void biquad_setQ(biquad* b, double Q);
extern void biquad_setPeakGain(biquad* b, double peakGain);

// get filter parameters
extern double biquad_getFs(const biquad* b);
extern bqFilterType biquad_getType(const biquad* b);
extern double biquad_getFc(const biquad* b);
extern double biquad_getQ(const biquad* b);
extern double biquad_getPeakGain(const biquad* b);

// clear buffers
extern void biquad_clear(biquad* b);

// return the filter name
const char* biquad_getFilterName(const biquad* b);

#ifdef __cplusplus
}
#endif

#endif // BIQUADS_H_
