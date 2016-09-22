#ifndef BIQUADS_HPP_
#define BIQUADS_HPP_

#include "biquads.h"

class biquadFilter {

public:

    // constructor
    biquadFilter(
        bqFilterType type = BQ_LOWPASS,
        int fs = 44100,
        double fc = 1000,
        double Q = 0.7071,
        double peakGain = 0
    ) {
        b = biquad_create(type, fs, fc, Q, peakGain);
    }
    
    // destructor
    ~biquadFilter(void) { biquad_destroy(b); }
    
    // calculate magnitude
    double magnitude(double f) { return biquad_magnitude(b, f); }
    
    // process sample
    bqFloat process(bqFloat x) { return biquad_process(b, x); }
    
    // set filter parameter - fs
    void setFs(double fs) { return biquad_setFs(b, fs); }
    
    // set filter parameter - type
    void setType(bqFilterType type) { biquad_setType(b, type); }
    
    // set filter parameter - fc
    void setFc(double fc) { biquad_setFc(b, fc); }
    
    // set filter parameter - Q
    void setQ(double Q) { biquad_setQ(b, Q); }
    
    // set filter parameter - peakGain
    void setPeakGain(double peakGain) { biquad_setPeakGain(b, peakGain); }
    
    // clear buffers
    void clear(void) { biquad_clear(b); }

private:
    
    biquad *b;

};

#endif // BIQUADS_HPP_
