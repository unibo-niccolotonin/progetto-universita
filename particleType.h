#ifndef PARTICLE_TYPE_H
#define PARTICLE_TYPE_H

#include <string>

class ParticleType
{
    public:
    virtual std::string getName() const;
    virtual double getCharge() const;
    virtual double getMass() const;
    virtual double getWidth() const;
    
    virtual void print() const;

    ParticleType(std::string name, double charge, double mass);
    virtual ~ParticleType() = default;

    protected:
    std::string fName;
    double const fCharge;
    double const fMass;

};

#endif