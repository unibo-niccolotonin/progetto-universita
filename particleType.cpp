#include "particleType.h"
#include <iostream>
#include <string>

ParticleType::ParticleType(std::string name, double charge, double mass)
: fName(name), fCharge(charge), fMass(mass)
{}

std::string ParticleType::getName() const
{
    return fName;
}

double ParticleType::getCharge() const
{
    return fCharge;
}

double ParticleType::getMass() const
{
    return fMass;
}

void ParticleType::print() const
{
    std::cout << "Name: " << fName << "\tCharge: " << fCharge << "\tMass: " << fMass << "\n";
}

double ParticleType::getWidth() const
{
    return 0;
}