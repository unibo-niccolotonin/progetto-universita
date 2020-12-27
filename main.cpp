#include "particle.h"
#include "generate.h"
#include <vector>
#include <iostream>

int main()
{


    Particle::addParticleType("Pi+", 1, 0.13957);
    Particle::addParticleType("Pi-", -1, 0.13957);

    Particle::addParticleType("K+", 1, 0.49367);
    Particle::addParticleType("K-", -1, 0.49367);

    Particle::addParticleType("p+", 1, 0.93827);
    Particle::addParticleType("p-", 1, 0.93827);

    Particle::addParticleType("K*", 0, 0.89166, 0.050);

    int numEvents;
    int numParticles;

    std::cout << "Insert numEvents\n";

    std::cin >> numEvents;

    std::cout << "Insert numParticles\n";

    std::cin >> numParticles;

    generate(numEvents, numParticles);

    std::cout << "Done!\n";

    return 0;
}