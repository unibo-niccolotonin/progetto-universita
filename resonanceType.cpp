#include "resonanceType.h"
#include <string>
#include <iostream>

ResonanceType::ResonanceType(std::string name, double charge, double mass, double width) :
ParticleType(name, charge, mass), fWidth(width)
{}

double ResonanceType::getWidth() const
{
    return fWidth;
}

void ResonanceType::print() const
{
    std::cout << "Name: " << fName << "\tCharge: " << fCharge << "\tMass: " << fMass << "\tWidth: " << fWidth << "\n";
}