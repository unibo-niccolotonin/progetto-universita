#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TStyle.h"
#include <iostream>


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

    //Check if the numbers of particle types agree with the a priori population ratio

    std::cout << "Bin errors:\n";

    for (int i = 1; i < 11; i++)
    {
        std::cout << i <<": " << typeHisto->GetBinError(i) << "\n";
    }

    // Angular functions fit

    phiHisto->Fit("pol0");
    phiHisto->SetFillColor(kBlue);
    phiHisto->GetXaxis()->SetTitle("rad");

    thetaHisto->Fit("pol0");
    thetaHisto->SetFillColor(kBlue);
    thetaHisto->GetXaxis()->SetTitle("rad");
    // Impulse fit

    impulseHisto->Fit("expo");
    impulseHisto->SetFillColor(kBlue);
    impulseHisto->GetXaxis()->SetTitle("GeV/c^2");

    std::cout << "impulse histogram chi squared: " << impulseHisto->GetFunction("expo")->GetChisquare() << "\n";
    std::cout << "Impulse histogram mean: " << impulseHisto->GetFunction("expo")->GetParameter(0) << "\n";

    // Decayed fit

    invMassDecayHisto->Fit("gaus", "", "", 0.7, 1.1);
    invMassDecayHisto->SetFillColor(kBlue);
    invMassDecayHisto->GetXaxis()->SetTitle("GeV/c^2");

    TH1D* firstComparisonHisto = (TH1D*)invMassDifChargeKaonPionHisto->Clone();
    firstComparisonHisto->SetName("firstComparisonHisto");
    firstComparisonHisto->SetTitle("K* invariant mass Distribution (Kaon-Pion)");
    firstComparisonHisto->GetXaxis()->SetTitle("GeV/c^2");
    firstComparisonHisto->SetFillColor(kBlue);


	firstComparisonHisto->Add(invMassSameChargeKaonPionHisto, -1);
    


    TH1D* secondComparisonHisto = (TH1D*)invMassDifChargeHisto->Clone();
    secondComparisonHisto->SetName("secondComparisonHisto");
    secondComparisonHisto->SetTitle("K* invariant mass Distribution");
    secondComparisonHisto->GetXaxis()->SetTitle("GeV/c^2");
    secondComparisonHisto->SetFillColor(kBlue);
    
	secondComparisonHisto->Add( invMassSameChargeHisto, -1);


    auto* f1 = new TF1("f1", "[0] + TMath::gaus([1], [2])", 0, 5);

	firstComparisonHisto->Fit("gaus","", "", 0.7, 1.1);
    secondComparisonHisto->Fit("gaus","", "", 0.7, 1);


	analyze_file->Write();

    analyze_file->Close();

    file->Close();
}
