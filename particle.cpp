#include "particle.h"
#include "resonanceType.h"
#include <memory>
#include <string>
#include <iostream>
#include <map>

#include <cmath>
#include <cstdlib>

using PPointer = std::unique_ptr<ParticleType>;

std::map<std::string, PPointer> Particle::particleTypeTable;

Particle::Particle(std::string type, double Px, double Py, double Pz)
{
    if (!particleTypeTable.count(type))
    {
        std::cout << "The type you've entered isn't right you dummy!";

        throw 0; // da sostituire
    }

    fType = type;
    fPx = Px;
    fPy = Py;
    fPz = Pz;
}

int Particle::fNParticleType = 0;

void Particle::addParticleType(std::string name, double charge, double mass)
{
    
    if (particleTypeTable.count(name))
    {
        std::cout << "You've already added a particle called " << name << " dummy!";
        throw 0; // da cambiare
    }

    ++fNParticleType;

    if (fNParticleType == fMaxNumParticleType)
    {
        std::cout << "You've added more than " << fMaxNumParticleType << " dummy!";
        throw 0; // da cambiare
    }

    //auto a = std::make_unique<ParticleType>(name, charge, mass); 
    particleTypeTable.emplace(name, std::make_unique<ParticleType>(name, charge, mass));
}

void Particle::addParticleType(std::string name, double charge, double mass, double width)
{
    
    if (particleTypeTable.count(name))
    {
        std::cout << "You've already added a particle called " << name << " dummy!";
        throw 0; // da cambiare
    }

    ++fNParticleType;

    if (fNParticleType == fMaxNumParticleType)
    {
        std::cout << "You've added more than " << fMaxNumParticleType << " dummy!";
        throw 0; // da cambiare
    }

    //auto a = std::make_unique<ParticleType>(name, charge, mass); 
    particleTypeTable.emplace(name, std::make_unique<ResonanceType>(name, charge, mass, width));
}

void Particle::printTypes()
{
    for (auto it = Particle::particleTypeTable.cbegin(); it != Particle::particleTypeTable.cend(); ++it)
    {
        it->second->print();
    }
}

void Particle::print()
{
    std::cout << "Name: " << fType << "\tPx: " << fPx << "\tPy: " << fPy << "\tPz: " << fPz << "\n"; 
}

double Particle::getMass() const
{
    return Particle::particleTypeTable.at(fType)->getMass();
}

double Particle::getCharge() const
{
    return Particle::particleTypeTable.at(fType)->getCharge();
}

double Particle::getEnergy() const
{
    return std::sqrt(getMass() * getMass() + 
    fPx * fPx + 
    fPy * fPy + 
    fPz * fPz);
}

double Particle::invMass(Particle const& other) const
{
    return std::sqrt((getEnergy() + other.getEnergy()) * (getEnergy() + other.getEnergy()) -
    (fPx + other.fPx) * (fPx + other.fPx) -
    (fPy + other.fPy) * (fPy + other.fPy) -
    (fPz + other.fPz) * (fPz + other.fPz));
}

std::string Particle::getType() const
{
    return fType;
}

double Particle::getPx() const
{
    return fPx;
}

double Particle::getPy() const
{
    return fPy;
}

double Particle::getPz() const
{
    return fPz;
}


void Particle::setType(std::string type)
{
    if(particleTypeTable.count(type))
    {
        fType = type;
    } else 
    {
        std::cout << "The type You specified (" << type << ") doesn't exist you dummy!\n";
        throw 2;
    }

}

void Particle::setP(double Px, double Py, double Pz)
{
    fPx = Px;
    fPy = Py;
    fPz = Pz;
}

int Particle::Decay2body(Particle &dau1,Particle &dau2) const {
  if(getMass() == 0.0){
    std::cout << "Decayment cannot be preformed if mass is zero\n";
    return 1;
  }
  
  double massMot = getMass();
  double massDau1 = dau1.getMass();
  double massDau2 = dau2.getMass();

  if(particleTypeTable.count(fType)){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += particleTypeTable.at(fType)->getWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    std::cout << "Decayment cannot be preformed because mass is too low in this channel\n";
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.setP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.setP(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(fPx*fPx + fPy*fPy + fPz*fPz + massMot*massMot);

  double bx = fPx/energy;
  double by = fPy/energy;
  double bz = fPz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz)
{

  double energy = getEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*fPx + by*fPy + bz*fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  fPx += gamma2*bp*bx + gamma*bx*energy;
  fPy += gamma2*bp*by + gamma*by*energy;
  fPz += gamma2*bp*bz + gamma*bz*energy;
}
