
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include <string>
#include <array>
#include <vector>
#include <map>
#include <dirent.h>
#include <iostream>
#include "TPaveText.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include <bits/stdc++.h>
#include <iostream>
#include<chrono>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include "ROOT/RDF/Utils.hxx"

TGraphErrors* sigmaDCAvsPt(ROOT::RDF::RResultPtr<TH2D> hDCAvsPt, std::string resultName, std::string axis, std::string dir){


    std::unique_ptr<TFile> projections(TFile::Open(Form("projections_%s.root",dir.c_str()), "RECREATE"));
    double sigma[50];
    double sigmaErrors[50];
    double pt[50];

    for (int i = 2; i < 51; i++){

        TH1D *h = hDCAvsPt->ProjectionY(Form("bin%d",i), i,i);
        TF1 fitF = TF1(Form("fitF_%d",i), "gaus", -0.02, 0.02);
        h->Fit(Form("fitF_%d",i), "R+ Q");
        TF1 fitF2 = TF1(Form("fitF2_%d",i), "gaus", -fitF.GetParameter(2) * 2, fitF.GetParameter(2) * 2);
        h->Fit(Form("fitF2_%d",i), "R+ Q");   /// refitting cutting the tails at 2*sigma
        // std::cout<<"Diff sigma: "<<fitF.GetParameter(2)<< "   "<<fitF2.GetParameter(2)<<std::endl;
        sigma[i-1] = fitF2.GetParameter(2)*10000;
        sigmaErrors[i-1] = fitF2.GetParError(2)*10000;
        pt[i-1] = (i-1)*0.2;
        h->Write();
    }
    // projections->Close();

    TGraphErrors* DCAresolution = new TGraphErrors(50, pt, sigma, 0, sigmaErrors) ;

    // TCanvas c;
    // DCAresolution->Draw("AP");
    std::string title = ("; p_{T} (GeV/c); #sigma_{DCA} (#mum)");
    DCAresolution->SetTitle(title.c_str());
    DCAresolution->SetMarkerColor(8);
    DCAresolution->SetMarkerStyle(20);
    DCAresolution->SetName(resultName.c_str());
    DCAresolution->GetXaxis()->SetRangeUser(0,10);
    DCAresolution->GetYaxis()->SetRangeUser(0,90);
    // std::unique_ptr<TGraphErrors> DCAresolution(50, pt, sigma, 0, sigmaErrors);
    // c.SaveAs(Form("%s.png", resultName.c_str()));

    // DCAresolution->Write();
    return DCAresolution;
};

template<typename T>
void fillHistograms(T table, TFile* myFile, std::string dir ){

    /// histograms 1D: masses
    auto hMassOmega = table.Histo1D({"hMassOmega", ";Invariant mass (GeV/#it{c}^2);Counts", 125, 1.600, 1.800},"fMassOmega");
    auto hMassV0 = table.Histo1D({"hMassV0", ";Invariant mass (GeV/#it{c}^2);Counts", 125, 1.090, 1.140},"fMassV0");
    auto hMassXi = table.Histo1D({"hMassXi", ";Invariant mass (GeV/#it{c}^2);Counts", 125, 1.310, 1.335},"fMassXi");

    /// fitting masses:
    //omega
    TF1 fitMassOmega = TF1("fitMassOmega", "gaus", 1.668, 1.677);
    fitMassOmega.SetLineWidth(0);
    TF1 fitBgOmega = TF1("fitBgOmega", "[0]*exp(-[1]*x)", 1.620, 1.800);
    fitBgOmega.SetLineWidth(0);
    fitBgOmega.SetParNames("a", "b");
    fitBgOmega.SetParameters(1, 1);
    hMassOmega->Fit("fitMassOmega", "R+ Q");
    hMassOmega->Fit("fitBgOmega", "R+ Q");    
    double parOmega[6];
    fitMassOmega.GetParameters(&parOmega[0]);
    fitBgOmega.GetParameters(&parOmega[3]);
    TF1 fitOmega = TF1("fitOmega", "fitMassOmega+fitBgOmega", 1.620, 1.800);
    fitOmega.SetParameters(parOmega);
    hMassOmega->Fit("fitOmega", "R+ Q");
    //v0
    TF1 fitMassV0 = TF1("fitMassV0", "gaus", 1.105, 1.125);
    hMassV0->Fit("fitMassV0", "R Q");
    //xi
    TF1 fitMassXi = TF1("fitMassXi", "gaus", 1.296, 1.346);
    hMassXi->Fit("fitMassXi", "R Q");

    /// histograms 1D: DCA
    auto hCascDCAxy = table.Histo1D({"hCascDCAxy", ";Cascade DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
    auto hCascDCAz = table.Histo1D({"hCascDCAz", ";Cascade DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");
    auto hProtonDCAxy = table.Histo1D({"hProtonDCAxy", ";Proton DCA_{xy} (cm);Counts", 200, -5., 5.},"fProtonDCAxy");
    auto hProtonDCAz = table.Histo1D({"hProtonDCAz", ";Proton DCA_{z} (cm);Counts", 200, -5., 5.},"fProtonDCAz");
    auto hPionDCAxy = table.Histo1D({"hPionDCAxy", ";Pion DCA_{xy} (cm);Counts", 200, -10., 10.},"fPionDCAxy");
    auto hPionDCAz = table.Histo1D({"hPionDCAz", ";Pion DCA_{z} (cm);Counts", 200, -10., 10.},"fPionDCAz");
    auto hBachDCAxy = table.Histo1D({"hBachDCAxy", ";Bachelor DCA_{xy} (cm);Counts", 200, -10., 10.},"fBachDCAxy");
    auto hBachDCAz = table.Histo1D({"hBachDCAz", ";Bachelor DCA_{z} (cm);Counts", 200, -10., 10.},"fBachDCAz");


    /// histograms 2D: mass vs DCA
    auto hOmegaMassVsDCAxy = table.Histo2D({"hOmegaMassVsDCAxy", ";#Omega DCA_{xy} (cm);Invariant mass (GeV/#it{c}^2)", 100, -0.02, 0.02, 100, 1.600, 1.800}, "fCascDCAxy", "fMassOmega");
    auto hOmegaMassVsDCAz = table.Histo2D({"hOmegaMassVsDCAz", ";#Omega DCA_{z} (cm);Invariant mass (GeV/#it{c}^2)", 100, -0.02, 0.02, 100, 1.600, 1.800}, "fCascDCAz", "fMassOmega");
    auto hXiMassVsDCAxy = table.Histo2D({"hXiMassVsDCAxy", ";#Xi DCA_{xy} (cm);Invariant mass (GeV/#it{c}^2)", 100, -0.02, 0.02, 100, 1.300, 1.400}, "fCascDCAxy", "fMassXi");
    auto hXiMassVsDCAz = table.Histo2D({"hXiMassVsDCAz", ";#Xi DCA_{z} (cm);Invariant mass (GeV/#it{c}^2)", 100, -0.02, 0.02, 100, 1.300, 1.400}, "fCascDCAz", "fMassXi");

    /// histograms 2D: DCA vs pt (no cuts)
    auto hDCAxyVsPt = table.Histo2D({"hDCAxyVsPt", ";p_{T} (GeV/#it{c});DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
    auto hDCAzVsPt = table.Histo2D({"hDCAzVsPt", ";p_{T} (GeV/#it{c});DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");
  
    TGraphErrors * hDCAxyResolution = sigmaDCAvsPt(hDCAxyVsPt, "hDCAxyResolution", "xy", dir.substr(0, dir.size()-1));
    TGraphErrors * hDCAzResolution = sigmaDCAvsPt(hDCAzVsPt, "hDCAzResolution","z", dir.substr(0, dir.size()-1));
    
    /// histograms 2D: DCA vs pt (mass cut: 2 sigma)
    // xi
    double meanMassXi = fitMassXi.GetParameter(1);
    double sigmaMassXi = fitMassXi.GetParameter(2);
    auto casctableCutXiMass = table.Filter([&](double fMassXi){if (fMassXi > meanMassXi - 2. * sigmaMassXi && fMassXi < meanMassXi + 2. * sigmaMassXi) return true; else return false;},{"fMassXi"});
    auto hXiCutsDCAxyVsPt = casctableCutXiMass.Histo2D({"hXiCutsDCAxyVsPt", ";p_{T} (GeV/#it{c});#Xi DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
    auto hXiCutsDCAzVsPt = casctableCutXiMass.Histo2D({"hXiCutsDCAzVsPt", ";p_{T} (GeV/#it{c});#Xi DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");
    //omega
    double meanMassOmega = fitMassOmega.GetParameter(1);
    double sigmaMassOmega = fitMassOmega.GetParameter(2);
    auto casctableCutOmegaMass = table.Filter([&](double fMassOmega){if (fMassOmega > meanMassOmega - 2. * sigmaMassOmega && fMassOmega < meanMassOmega + 2. * sigmaMassOmega) return true; else return false;},{"fMassOmega"});
    auto hOmegaCutsDCAxyVsPt = casctableCutOmegaMass.Histo2D({"hOmegaCutsDCAxyVsPt", ";p_{T} (GeV/#it{c});#Omega DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
    auto hOmegaCutsDCAzVsPt = casctableCutOmegaMass.Histo2D({"hOmegaCutsDCAzVsPt", ";p_{T} (GeV/#it{c});#Omega DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");


    // histograms 2D: mass vs pt
    auto hXiMassVsPt = table.Histo2D({"hXiMassVsPt", ";p_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^2)", 50, 0.,10., 125, 1.310, 1.335}, "fCascPt", "fMassXi");
    auto hOmegaMassVsPt = table.Histo2D({"hOmegaMassVsPt", ";p_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^2)", 50, 0.,10., 125, 1.600, 1.800}, "fCascPt", "fMassOmega");
    
    // applying PID cuts
    
    myFile->cd(dir.c_str());
    hMassOmega->Write();
    hMassV0->Write();
    hMassXi->Write();

    hCascDCAxy->Write();
    hCascDCAz->Write();
    hProtonDCAxy->Write();
    hProtonDCAz->Write();
    hPionDCAxy->Write();
    hPionDCAz->Write();
    hBachDCAxy->Write();
    hBachDCAz->Write();

    hOmegaMassVsDCAxy->Write();
    hOmegaMassVsDCAz->Write();
    hXiMassVsDCAxy->Write();
    hXiMassVsDCAz->Write();

    hDCAxyVsPt->Write();
    hDCAzVsPt->Write();

    hDCAxyResolution->Write();
    hDCAzResolution->Write();

    hXiCutsDCAxyVsPt->Write();
    hXiCutsDCAzVsPt->Write();
    hOmegaCutsDCAxyVsPt->Write();
    hOmegaCutsDCAzVsPt->Write();

    hXiMassVsPt->Write();
    hOmegaMassVsPt->Write();
    
}





void readAO2D(){

    int minClusTPC = 70;
    float minNsigmaTPC = -4.;
    float maxNsigmaTPC = 4.;


    ROOT::RDataFrame casctable("DF_2316989016393920/O2npcasctable", "AO2D_apass2.root");

    auto omegatable = casctable.Filter("abs(fMassXi-1.32171) > 0.005"); //(massXi - 1.32171) > 0.005
    auto casctablePIDcuts = casctable.Filter("fProtonNClusTPC>70 and fPionNClusTPC>70 and fProtonTPCNSigma > -4 and fProtonTPCNSigma<4 and fPionTPCNSigma>-4 and fPionTPCNSigma<4 and fBachKaonTPCNSigma>-4 and fBachKaonTPCNSigma<4 and fBachPionTPCNSigma>-4 and fBachPionTPCNSigma<4");
    auto omegatablePIDcuts = omegatable.Filter("fProtonNClusTPC>70 and fPionNClusTPC>70 and fProtonTPCNSigma > -4 and fProtonTPCNSigma<4 and fPionTPCNSigma>-4 and fPionTPCNSigma<4 and fBachKaonTPCNSigma>-4 and fBachKaonTPCNSigma<4");
    
    
    // auto casctablePIDcuts = casctable.Filter("fProtonNClusTPC>70 and fPionNClusTPC>70  and fProtonTPCNSigma > -4 and fProtonTPCNSigma4  ");
    // auto omegatablePIDcuts = omegatable.Filter("fProtonNClusTPC>70 and fPionNClusTPC>70 and fProtonTPCNSigma > -4 and fProtonTPCNSigma>4 and fPionTPCNSigma>-4 and fPionTPCNSigma<4 and fBachKaonTPCNSigma>-4 and fBachKaonTPCNSigma<4");




    // std::unique_ptr<TFile> myFile(TFile::Open("histos.root", "RECREATE"));
    TFile* myFile(TFile::Open("histos.root", "RECREATE"));
    myFile->mkdir("Cascade/");
    myFile->mkdir("Omega/");
    myFile->mkdir("Cascade+PIDcuts/");
    myFile->mkdir("Omega+PIDcuts/");
    


    fillHistograms(casctable, myFile, "Cascade/");
    fillHistograms(omegatable, myFile, "Omega/");
    fillHistograms(casctablePIDcuts, myFile, "Cascade+PIDcuts/");
    fillHistograms(omegatablePIDcuts, myFile, "Omega+PIDcuts/");

    myFile->Close();

}