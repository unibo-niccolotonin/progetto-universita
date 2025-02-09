#include "generate.h"
#include "particle.h"
#include "TRandom.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include <vector>
#include <iostream>


void generate(int numEvents, int numParticles)
{

    gRandom->SetSeed();

    std::vector<Particle> particleVector(numParticles + numParticles/2, Particle("K*", 0, 0, 0));

    int n_bins = 80;

    auto* ALICEFile = new TFile("ALICE.root", "RECREATE");

    auto* phiHisto = new TH1D("phiHisto", "phiHisto", n_bins, 0, 2 * TMath::Pi());
    auto* thetaHisto = new TH1D("thetaHisto", "thetaHisto", n_bins, 0, TMath::Pi());
    auto* typeHisto = new TH1I("typeHisto", "typehisto", 10, 0, 9);
    auto* impulseHisto = new TH1D("impulseHisto", "impulseHisto", n_bins, 0, 5);
    auto* energyHisto = new TH1D("energyHisto", "energyHisto", n_bins, 0, 5);
    
    auto* invMassSameChargeHisto = new TH1D("invMassSameChargeHisto", "invMassSameChargeHisto", n_bins, 0, 5);
    auto* invMassDifChargeHisto = new TH1D("invMassDifChargeHisto", "invMassDifChargeHisto", n_bins, 0, 5);
    auto* invMassSameChargeKaonPionHisto = new TH1D("invMassSameChargeKaonPionHisto", "invMassSameChargeKaonPionHisto", n_bins, 0, 5);
    auto* invMassDifChargeKaonPionHisto = new TH1D("invMassDifChargeKaonPionHisto", "invMassDifChargeKaonPionHisto", n_bins, 0, 5);
    auto* invMassDecayHisto = new TH1D("invMassDecayHisto", "invMassDecayHisto", n_bins, 0, 5);

    {
        int index = 1;
        for (auto it = Particle::particleTypeTable.cbegin(); it != Particle::particleTypeTable.cend(); ++it)
        {
        typeHisto->GetXaxis()->SetBinLabel(index ,it->second->getName().c_str());
        ++index;
        }
    }

    for(int i = 0; i < numEvents; ++i)
    {
        int decayedIndex = 0;
        for (int j = 0; j < numParticles; ++j)
        {
            double random = gRandom->Uniform(0, 1);
            double phi = gRandom->Uniform(0, 2 * TMath::Pi());
            double theta = gRandom->Uniform(0, TMath::Pi());
            double impulse = gRandom->Exp(1);

            thetaHisto->Fill(theta);
            phiHisto->Fill(phi);
            impulseHisto->Fill(impulse);

            //Define momentum
            double Px = impulse * TMath::Sin(theta) * TMath::Cos(phi);
            double Py = impulse * TMath::Sin(theta) * TMath::Sin(phi);
            double Pz = impulse * TMath::Cos(theta);

            particleVector[j].setP(Px, Py, Pz);

            //Set ParticleType and momentum
            if (random < 0.4)
            {
                particleVector[j].setType("K+");
            } else if (random < 0.8)
            {
                particleVector[j].setType("K-");
            } else if (random < 0.85)
            {
                particleVector[j].setType("Pi+");
            } else if (random < 0.9)
            {
                particleVector[j].setType("Pi-");
            } else if (random < 0.945)
            {
                particleVector[j].setType("p+");
            } else if (random < 0.99)
            {
                particleVector[j].setType("p-");
            } else
            {
                particleVector[j].setType("K*");
                
                double rand = gRandom->Uniform(0, 1);

                if (rand <= 0.5)
                {
                    Particle &dau1 = particleVector[numParticles + decayedIndex];
                    Particle &dau2 = particleVector[numParticles + decayedIndex + 1];

                    dau1.setType("Pi+");
                    dau2.setType("K-");

                    particleVector[j].Decay2body(dau1, dau2);

                    invMassDecayHisto->Fill(dau1.invMass(dau2));

                    decayedIndex += 2;
                } else 
                {
                    Particle &dau1 = particleVector[numParticles + decayedIndex];
                    Particle &dau2 = particleVector[numParticles + decayedIndex + 1];

                    dau1.setType("Pi-");
                    dau2.setType("K+");

                    particleVector[j].Decay2body(dau1, dau2);

                    invMassDecayHisto->Fill(dau1.invMass(dau2));

                    decayedIndex += 2;
                }
                
            }

            energyHisto->Fill(particleVector[j].getEnergy());
            typeHisto->Fill(particleVector[j].getType().c_str(), 1);
        }

        //Fill invariant mass histograms
        for (int k = 0; k < (numParticles + decayedIndex); ++k)
        {

            if (particleVector[k].getCharge() == 0) continue;

            for (int l = (k+1); l < (numParticles + decayedIndex); ++l)
            {
                int chargeSignAgreement = particleVector[k].getCharge() * particleVector[l].getCharge();

                switch (chargeSignAgreement)
                {
                case 1:
                    invMassSameChargeHisto->Fill(particleVector[k].invMass(particleVector[l]));
                    if (particleVector[k].getType()[0] == 'K' && particleVector[l].getType()[0] == 'P' ||
                        particleVector[k].getType()[0] == 'P' && particleVector[l].getType()[0] == 'K')
                        invMassSameChargeKaonPionHisto->Fill(particleVector[k].invMass(particleVector[l]));
                    break;
                
                case -1:
                    invMassDifChargeHisto->Fill(particleVector[k].invMass(particleVector[l]));
                    if (particleVector[k].getType()[0] == 'K' && particleVector[l].getType()[0] == 'P' ||
                        particleVector[k].getType()[0] == 'P' && particleVector[l].getType()[0] == 'K')
                        invMassDifChargeKaonPionHisto->Fill(particleVector[k].invMass(particleVector[l]));
                    break;

                default:
                    break;
                }
            }
        }
        


        decayedIndex = 0;
        std::cout << i << "\n";
    }

    ALICEFile->Write();
    
    ALICEFile->Close();

}
