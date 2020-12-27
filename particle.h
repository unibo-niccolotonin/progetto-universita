#ifndef PARTICLE_H
#define PARTICLE_H

#include "particleType.h"
#include <string>
#include <map>
#include <memory>

class Particle
{
    using PPointer = std::unique_ptr<ParticleType>;


    public:
    Particle(std::string type, double Px = 0, double Py = 0, double Pz = 0);

    static int fNParticleType;
    static void addParticleType(std::string name, double charge, double mass);
    static void addParticleType(std::string name, double charge, double mass, double width);
    static void printTypes();
    static std::map<std::string, PPointer> particleTypeTable;     //trovare modo per rendere privato

    double getMass() const;
    double getCharge() const;
    double getEnergy() const;
    double invMass(Particle const& other) const;

    std::string getType() const;
    double getPx() const;
    double getPy() const;
    double getPz() const;

    void setType(std::string type);
    void setP(double Px, double Py, double Pz);

    int Decay2body(Particle &dau1,Particle &dau2) const;

    void print();

    private:

    void Boost(double bx, double by, double bz);

    std::string fType;
    double fPx;
    double fPy;
    double fPz;

    static int const fMaxNumParticleType = 10;
};

#endif