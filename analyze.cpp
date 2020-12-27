#include "TFile.h"
#include "TH1.h"
#include <iostream>

double chi_squared_bin(double totalEntries, double observed, double probability)
{
    double expected = totalEntries * probability;
    double sigmaSquared = expected * (1 - probability);

    return ((observed - expected) * (observed - expected) / sigmaSquared);

}

void analyze()
{
    auto* file = new TFile("ALICE.root");

    TH1D* phiHisto = (TH1D*)file->Get("phiHisto");
    TH1D* thetaHisto = (TH1D*)file->Get("thetaHisto");
    TH1I* typeHisto = (TH1I*)file->Get("typeHisto");
    TH1D* impulseHisto = (TH1D*)file->Get("impulseHisto");
    TH1D* energyHisto = (TH1D*)file->Get("energyHisto");

    TH1D* invMassSameChargeHisto = (TH1D*)file->Get("invMassSameChargeHisto");
    TH1D* invMassDifChargeHisto = (TH1D*)file->Get("invMassDifChargeHisto");
    TH1D* invMassSameChargeKaonPionHisto = (TH1D*)file->Get("invMassSameChargeKaonPionHisto");
    TH1D* invMassDifChargeKaonPionHisto = (TH1D*)file->Get("invMassDifChargeKaonPionHisto");
    TH1D* invMassDecayHisto = (TH1D*)file->Get("invMassDecayHisto");

    std::cout << "Number of entries in the histograms\n";

    std::cout << "phiHisto:" << phiHisto->GetEntries() << "\n";
    std::cout << "thetaHisto:" << thetaHisto->GetEntries() << "\n";
    std::cout << "typeHisto:" << typeHisto->GetEntries() << "\n";
    std::cout << "impulseHisto:" << impulseHisto->GetEntries() << "\n";
    std::cout << "energyHisto:" << energyHisto->GetEntries() << "\n";

    std::cout << "invMassSameChargeHisto:" << invMassSameChargeHisto->GetEntries() << "\n";
    std::cout << "invMassDifChargeHisto:" << invMassDifChargeHisto->GetEntries() << "\n";
    std::cout << "invMassSameChargeKaonPionHisto:" << invMassSameChargeKaonPionHisto->GetEntries() << "\n";
    std::cout << "invMassDifChargeKaonPionHisto:" << invMassDifChargeKaonPionHisto->GetEntries() << "\n";
    std::cout << "invMassDecayHisto:" << invMassDecayHisto->GetEntries() << "\n";

    //Check if the numbers of particle types agree with the a priori population ratio

    auto totalTypeEntries = typeHisto->GetEntries();

    std::cout << "Chi squared of type bins:\n";

    std::cout << "Pi+: " << chi_squared_bin(totalTypeEntries, typeHisto->GetBinContent(typeHisto->GetXaxis()->FindBin("Pi+")), 0.4) << "\n";
    std::cout << "Pi-: " << chi_squared_bin(totalTypeEntries, typeHisto->GetBinContent(typeHisto->GetXaxis()->FindBin("Pi-")), 0.4) << "\n";

    std::cout << "K+: " << chi_squared_bin(totalTypeEntries, typeHisto->GetBinContent(typeHisto->GetXaxis()->FindBin("K+")), 0.05) << "\n";
    std::cout << "K-: " << chi_squared_bin(totalTypeEntries, typeHisto->GetBinContent(typeHisto->GetXaxis()->FindBin("K-")), 0.05) << "\n";

    std::cout << "p+: " << chi_squared_bin(totalTypeEntries, typeHisto->GetBinContent(typeHisto->GetXaxis()->FindBin("p+")), 0.045) << "\n";
    std::cout << "p-: " << chi_squared_bin(totalTypeEntries, typeHisto->GetBinContent(typeHisto->GetXaxis()->FindBin("p-")), 0.045) << "\n";

    std::cout << "K*: " << chi_squared_bin(totalTypeEntries, typeHisto->GetBinContent(typeHisto->GetXaxis()->FindBin("K*")), 0.01) << "\n";

    //
    // Implementare la stampatura degli errori
    //

    // Angular functions fit

    phiHisto->Fit("pol0");
    std::cout << "phi histogram chi squared: " << phiHisto->GetFunction("pol0")->GetChisquare();

    thetaHisto->Fit("pol0");
    std::cout << "theta histogram chi squared: " << thetaHisto->GetFunction("pol0")->GetChisquare();

    // Impulse fit

    impulseHisto->Fit("expo");
    std::cout << "impulse histogram chi squared: " << impulseHisto->GetFunction("expo")->GetChisquare();
    std::cout << "Impulse histogram mean: " << impulseHisto->GetFunction("expo")->GetParameter(0);

    //
    // Implementare la stampatura degli errori
    //

    // Decayed fit

    auto* firstComparisonHisto = invMassSameChargeKaonPionHisto - invMassDifChargeKaonPionHisto;
    auto* firstDecayedFit = firstComparisonHisto->Fit(invMassDecayHisto);

    std::cout << "firstComparisonHisto chi squared: " << firstDecayedFit->GetChiSquared();

    auto* secondComparisonHisto = invMassSameChargeHisto - invMassDifChargeHisto;
    auto* secondDecayedFit = secondComparisonHisto->Fit(firstComparisonHisto);

    std::cout << "secondComparisonHisto chi squared: " << secondDecayedFit->GetChiSquared();

    auto* gaussianFit = firstComparisonHisto->Fit("gaus");
    std::cout << "K* Mass: " << gaussianFit->getParameter(0);
    std::cout << "K* Width: " << gaussianFit->getParameter(1);

    file->Close();
}