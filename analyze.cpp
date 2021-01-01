#include "TFile.h"
#include "TF1.h"
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
    gStyle->SetOptFit(1111);

    auto* file = new TFile("ALICE.root");
    auto* analyze_file = new TFile("ANALYSIS.root","RECREATE");

    TH1D* phiHisto = (TH1D*)file->Get("phiHisto")->Clone();
    TH1D* thetaHisto = (TH1D*)file->Get("thetaHisto")->Clone();
    TH1I* typeHisto = (TH1I*)file->Get("typeHisto")->Clone();
    TH1D* impulseHisto = (TH1D*)file->Get("impulseHisto")->Clone();
    TH1D* energyHisto = (TH1D*)file->Get("energyHisto")->Clone();

    TH1D* invMassSameChargeHisto = (TH1D*)file->Get("invMassSameChargeHisto")->Clone();
    TH1D* invMassDifChargeHisto = (TH1D*)file->Get("invMassDifChargeHisto")->Clone();
    TH1D* invMassSameChargeKaonPionHisto = (TH1D*)file->Get("invMassSameChargeKaonPionHisto")->Clone();
    TH1D* invMassDifChargeKaonPionHisto = (TH1D*)file->Get("invMassDifChargeKaonPionHisto")->Clone();
    TH1D* invMassDecayHisto = (TH1D*)file->Get("invMassDecayHisto")->Clone();

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

    using BinContent = typeHisto->GetBinContent;
    using Bin = typeHisto->GetXaxis()->FindBin;

    std::cout << "Pi+: " << chi_squared_bin(totalTypeEntries, BinContent(Bin("Pi+")), 0.4) << "\n";
    std::cout << "Pi-: " << chi_squared_bin(totalTypeEntries, BinContent(Bin("Pi-")), 0.4) << "\n";

    std::cout << "K+: " << chi_squared_bin(totalTypeEntries, BinContent(Bin("K+")), 0.05) << "\n";
    std::cout << "K-: " << chi_squared_bin(totalTypeEntries, BinContent(Bin("K-")), 0.05) << "\n";

    std::cout << "p+: " << chi_squared_bin(totalTypeEntries, BinContent(Bin("p+")), 0.045) << "\n";
    std::cout << "p-: " << chi_squared_bin(totalTypeEntries, BinContent(Bin("p-")), 0.045) << "\n";

    std::cout << "K*: " << chi_squared_bin(totalTypeEntries, BinContent(Bin("K*")), 0.01) << "\n";

    //
    // Implementare la stampatura degli errori
    //

    // Angular functions fit

    phiHisto->Fit("pol0");
    std::cout << "phi histogram chi squared: " << phiHisto->GetFunction("pol0")->GetChisquare() << "\n";

    thetaHisto->Fit("pol0");
    std::cout << "theta histogram chi squared: " << thetaHisto->GetFunction("pol0")->GetChisquare() << "\n";

    // Impulse fit

    impulseHisto->Fit("expo");
    std::cout << "impulse histogram chi squared: " << impulseHisto->GetFunction("expo")->GetChisquare() << "\n";
    std::cout << "Impulse histogram mean: " << impulseHisto->GetFunction("expo")->GetParameter(0) << "\n";

    //
    // Implementare la stampatura degli errori
    //

    // Decayed fit

    TH1D* firstComparisonHisto = (TH1D*)invMassDifChargeKaonPionHisto->Clone();
    firstComparisonHisto->SetName("firstComparisonHisto");
	firstComparisonHisto->Add(invMassSameChargeKaonPionHisto, -1); 

	double first_dec_fit = firstComparisonHisto->Chi2Test(invMassDecayHisto, "CHI2/NDF");

    std::cout << "firstComparisonHisto chi squared/NDF: " << first_dec_fit << "\n";

    TH1D* secondComparisonHisto = (TH1D*)invMassSameChargeHisto->Clone();
    secondComparisonHisto->SetName("secondComparisonHisto");
	secondComparisonHisto->Add( invMassDifChargeHisto, -1); 


	double second_dec_fit = secondComparisonHisto->Chi2Test(firstComparisonHisto,"CHI2/NDF");

    std::cout << "secondComparisonHisto chi squared/NDF: " << second_dec_fit << "\n";

	firstComparisonHisto->Fit("gaus");


    std::cout << "K* Mass: " << firstComparisonHisto->GetFunction("gaus")->GetParameter(0) << "\n";
    std::cout << "K* Width: " << firstComparisonHisto->GetFunction("gaus")->GetParameter(1) << "\n";


	analyze_file->Write();

    file->Close();
}
