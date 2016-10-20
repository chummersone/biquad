// A basic C++ wrapper for the biquad filter

#ifndef BIQUADS_HPP_
#define BIQUADS_HPP_

#include <memory>
#include <new>
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
    );
    
    // destructor
    ~biquadFilter(void) = default;
    
    // copy constructor
    biquadFilter(const biquadFilter& other);
    
    // move constructor
    biquadFilter(biquadFilter&& other);
    
    // copy assignment
    biquadFilter& operator=(biquadFilter other);
    
    // move assignment
    biquadFilter& operator=(biquadFilter&& other);
    
    // swap
    friend void swap(biquadFilter& first, biquadFilter& second);
    
    // calculate magnitude
    double magnitude(double f) const;
    
    // process sample
    bqFloat process(bqFloat x);
    
    // set/get filter parameter - fs
    void setFs(double fs);
    double getFs(void) const;
    
    // set/get filter parameter - type
    void setType(bqFilterType type);
    bqFilterType getType(void) const;
    
    // set/get filter parameter - fc
    void setFc(double fc);
    double getFc(void) const;
    
    // set/get filter parameter - Q
    void setQ(double Q);
    double getQ(void) const;
    
    // set/get filter parameter - peakGain
    void setPeakGain(double peakGain);
    double getPeakGain(void) const;
    
    // clear buffers
    void clear(void);

private:
    
    std::unique_ptr<biquad, void(*)(biquad*)> b;

};

// Implementation

// constructor
biquadFilter::biquadFilter(
    bqFilterType type,
    int fs,
    double fc,
    double Q,
    double peakGain
) : b(biquad_create(type, fs, fc, Q, peakGain), biquad_destroy) {
    if (!b) throw std::bad_alloc();
}

// copy constructor
biquadFilter::biquadFilter(const biquadFilter& other) :
b(biquad_clone(other.b.get()), biquad_destroy) {
    if (!b) throw std::bad_alloc();
}

// move constructor
biquadFilter::biquadFilter(biquadFilter&& other) :
biquadFilter() {
    swap(*this, other);
}

// copy assignment
biquadFilter& biquadFilter::operator=(biquadFilter other) {
    swap(*this, other);
    return *this;
}

// move assignment
biquadFilter& biquadFilter::operator=(biquadFilter&& other) {
    this->b = std::move(other.b);
    return *this;
}

// swap
void swap(biquadFilter& first, biquadFilter& second) {
    using std::swap;
    swap(first.b, second.b);
}

// calculate magnitude
double biquadFilter::magnitude(double f) const {
    return biquad_magnitude(b.get(), f);
}

// process sample
bqFloat biquadFilter::process(bqFloat x) {
    return biquad_process(b.get(), x);
}

// set/get filter parameter - fs
void biquadFilter::setFs(double fs) {
    return biquad_setFs(b.get(), fs);
}
double biquadFilter::getFs(void) const {
    return biquad_getFs(b.get());
}

// set/get filter parameter - type
void biquadFilter::setType(bqFilterType type) {
    biquad_setType(b.get(), type);
}
bqFilterType biquadFilter::getType(void) const {
return biquad_getType(b.get());
}

// set/get filter parameter - fc
void biquadFilter::setFc(double fc) {
    biquad_setFc(b.get(), fc);
}
double biquadFilter::getFc(void) const {
    return biquad_getFc(b.get());
}

// set/get filter parameter - Q
void biquadFilter::setQ(double Q) {
    biquad_setQ(b.get(), Q);
}
double biquadFilter::getQ(void) const {
    return biquad_getQ(b.get());
}

// set/get filter parameter - peakGain
void biquadFilter::setPeakGain(double peakGain) {
    biquad_setPeakGain(b.get(), peakGain);
}
double biquadFilter::getPeakGain(void) const {
    return biquad_getPeakGain(b.get());
}

// clear buffers
void biquadFilter::clear(void) {
    biquad_clear(b.get());
}

#endif // BIQUADS_HPP_
