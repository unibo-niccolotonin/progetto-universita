#ifndef RESONANCE_TYPE_H
#define RESONANCE_TYPE_H

#include "particleType.h"
#include <string>

class ResonanceType: public ParticleType
{
    public:

    ResonanceType(std::string name, double charge, double mass, double width);

    double getWidth() const override;

    void print() const override;

    private:
    double const fWidth;
};

#endif