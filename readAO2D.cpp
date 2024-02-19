
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
#include <iostream>
#include<chrono>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include "ROOT/RDF/Utils.hxx"
#include <cmath>

#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooPoisson.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooBifurGauss.h"


std::tuple<double, double, double, double> fitRoofit(TH1* h, float rangea, float rangeb, std::string name, std::string dir, std::string fittype, std::string particle, float pt, bool StoN = false, bool DCA = false, std::string DCAtype = ""){
    // using namespace RooFit;

    RooRealVar mass( "mass", "mass", rangea, rangeb); 
    RooRealVar mu("mu", "mean",  rangea, rangeb, "GeV/c^{2}");
    RooRealVar sigma("sigma", "width", 0.001, 0.01);

    
    RooRealVar a1("a1", "a1", 0., 5.);
    RooRealVar a2("a2", "a2", 0., 5.);
    RooRealVar n1("n1", "n1", 0, 10.);
    RooRealVar n2("n2", "n2", 0, 10.);

    RooRealVar c0("c0", "constant c0", -100., 0.);

    RooRealVar n_signal("n_signal", "n_signal",0, 30000);
    RooRealVar n_background("n_background", "n_background",0,  30000);

    RooCrystalBall cb("cb", "cb", mass, mu, sigma, a1, n1, a2, n2);

    RooGaussian gaus("gaus", "gaus", mass, mu, sigma);

    RooExponential exp("exp", "exp", mass, c0);

    RooDataHist dh("dh", "dh",mass, h);

    TCanvas* c = new TCanvas( "c", "c", 800, 600 );
    RooPlot* xFrame = mass.frame( RooFit::Title( "; Inv. Mass (GeV/c^{2}); " ) );

    if(DCA){
        xFrame->SetTitle(Form("; DCA_{%s} (cm); ", DCAtype.c_str()));
    }
    
    dh.plotOn(xFrame);
    
    if (particle == "Omega"){
        std::cout << "Fitting Omega" << std::endl;
        RooAddPdf modelOmega("modelOmega", "modelOmega", RooArgList(cb, exp), RooArgList(n_signal, n_background) );
        modelOmega.fitTo(dh, RooFit::Range(rangea, rangeb));
        modelOmega.plotOn(xFrame, RooFit::LineColor(kRed), RooFit::LineWidth(1));
    }

    else if(particle == "Xi" && StoN){
        std::cout << "Fitting Xi" << std::endl;
        RooAddPdf modelXi("modelXi", "modelXi", RooArgList(cb, exp), RooArgList(n_signal, n_background) );
        modelXi.fitTo(dh, RooFit::Range(rangea, rangeb));
        modelXi.plotOn(xFrame, RooFit::LineColor(kRed), RooFit::LineWidth(1));
    }

    if (fittype == "gaus" && particle != "Omega") {
        gaus.fitTo(dh, RooFit::Range(rangea, rangeb),RooFit::PrintLevel(-1));
    }
    else if (fittype == "cb" && particle != "Omega" && !StoN) {
        cb.fitTo(dh, RooFit::Range(rangea, rangeb), RooFit::PrintLevel(-1));
    }
    
    if (fittype == "gaus" && particle != "Omega") {
        gaus.plotOn(xFrame, RooFit::LineColor(kRed), RooFit::LineWidth(1)  );
    }
    else if (fittype == "cb" && particle != "Omega" && !StoN) {
        cb.plotOn(xFrame, RooFit::LineColor(kRed), RooFit::LineWidth(1)  );
    }
    
    xFrame->Draw();
    std::string filename = Form("plots/%s/%s", dir.c_str(), name.c_str());

    // std::cout<<"\n+++++++++++++++++++++++++++sigma "<<sigma.getVal()<<" "<<sigma.getError()<<" ++++++++++++++++"<<std::endl;

    TPaveText *pv = new TPaveText(0.6, 0.65, 0.8, 0.9, "NDC");
    pv->SetTextFont(132);
    pv->SetFillStyle(0);
    pv->SetBorderSize(0);
    pv->AddText(Form("#%s", particle.c_str()));
    ((TText*)pv->GetListOfLines()->Last())->SetTextSize(0.07);
    if (pt!=-1 && DCA){
        pv->AddText(Form("%0.3f < #it{p}_{T} < %0.3f GeV/#it{c}", pt-0.2, pt+0.2));
    }
    else if (pt!=-1){
        pv->AddText(Form("%0.3f < #it{p}_{T} < %0.3f GeV/#it{c}", pt-0.1, pt+0.1));
    }

    if (DCA){
        pv->AddText(Form("#mu = %.3f #pm %.3f #mum",  mu.getVal()*10000, mu.getError()*10000));
        pv->AddText(Form("#sigma = %.3f #pm %.3f #mum",  sigma.getVal()*10000,sigma.getError()*10000));
    }
    else{
        pv->AddText(Form("#mu = %.3f #pm %.3f GeV/#it{c}^{2}",  mu.getVal(), mu.getError()));
        pv->AddText(Form("#sigma = %.3f #pm %.3f GeV/#it{c}^{2}",  sigma.getVal(),sigma.getError()));

    }


    
    pv->Draw("same");
    c->SaveAs(filename.c_str());

    if (StoN){
        auto signal_counts = n_signal.getVal();
        auto signal_counts_errors = n_signal.getError();
        auto background_counts = n_background.getVal();
        auto background_counts_errors = n_background.getError();

        // signal within 3 sigma
        mass.setRange("signal", mu.getVal()-2*sigma.getVal(), mu.getVal()+2*sigma.getVal());
        auto signal_int = cb.createIntegral(RooArgSet(mass), RooArgSet(mass), "signal");
        auto signal_int_val_2s = signal_int->getVal()*signal_counts;
        auto signal_int_val_2s_error = signal_int_val_2s*signal_counts_errors/signal_counts;

        //background within 3 sigma
        mass.setRange("bkg", mu.getVal()-2*sigma.getVal(), mu.getVal()+2*sigma.getVal());
        auto bkg_int = exp.createIntegral(RooArgSet(mass), RooArgSet(mass), "bkg");
        auto bkg_int_val_2s = bkg_int->getVal()*background_counts;
        auto bkg_int_val_2s_error = bkg_int_val_2s*background_counts_errors/background_counts;

        // return std::tuple<double, double,double,double>{n_signal.getVal(),n_background.getVal(), n_signal.getError(), n_background.getError()};
        return std::tuple<double,double,double,double>{signal_int_val_2s,bkg_int_val_2s, signal_int_val_2s_error, bkg_int_val_2s_error};
    }
    if (DCA){
        return std::tuple<double, double,double,double>{mu.getVal(), sigma.getVal()*10000, mu.getError(), sigma.getError()*10000};
    }
    return std::tuple<double, double,double,double>{mu.getVal(), sigma.getVal(), mu.getError(), sigma.getError()};
}


TH1F*  stdMass(ROOT::RDF::RResultPtr<TH2D> hMassVsPt, std::string resultName, std::string dir, std::string particle, std::string title){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::unique_ptr<TFile> projections(TFile::Open(Form("projections/projectionsMass_%s.root",dir.c_str()), "RECREATE"));
    double pt[50];
    double sigma[50];
    double sigmaErrors[50];

    double mean[50];
    double meanError[50];

    for (int i = 0; i < 50; i++){   /// going from 1 GeV to 7 GeV

        TH1D *h = hMassVsPt->ProjectionY(Form("bin%d",i), i,i);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.5);

        if (particle == "Omega"){
            h->GetXaxis()->SetRangeUser(1.647, 1.697);

            std::tuple<double,double,double,double> mean_sigma_merr_serr = fitRoofit(h,1.647, 1.697, Form("stdMass/bin%d.pdf", i), dir , "cb", particle, 0.1+ i*0.2);
            
            pt[i] = 0.1 + i*0.2;
            sigma[i] = std::get<1>(mean_sigma_merr_serr);
            sigmaErrors[i] = std::get<3>(mean_sigma_merr_serr);
            mean[i] = std::get<0>(mean_sigma_merr_serr);
            meanError[i] = std::get<2>(mean_sigma_merr_serr);
        }

        else if (particle == "Xi"){
            std::tuple<double,double,double,double> mean_sigma_merr_serr = fitRoofit(h, 1.310, 1.335, Form("stdMass/bin%d.pdf", i), dir , "cb", particle, 0.1+ i*0.2);
            pt[i] = 0.1 + i*0.2;
            sigma[i] = std::get<1>(mean_sigma_merr_serr);
            sigmaErrors[i] = std::get<3>(mean_sigma_merr_serr);
            mean[i] = std::get<0>(mean_sigma_merr_serr);
            meanError[i] = std::get<2>(mean_sigma_merr_serr);
        }
        h->Write();
        h->SetStats(0);
    }
    
    projections->Close();

    TH1F* sdtMassVsPt = new TH1F("sdtMassVsPt", "", 50, 0, 10);
    sdtMassVsPt->SetMarkerStyle(21);
    sdtMassVsPt->SetTitle(title.c_str());
    sdtMassVsPt->SetName(resultName.c_str());
    sdtMassVsPt->SetMarkerStyle(21);
    sdtMassVsPt->GetXaxis()->SetRangeUser(0,10);
    sdtMassVsPt->GetYaxis()->SetRangeUser(0,10);
    for (int i=6; i<40; i++){    //between 1 and 8 GeV
        sdtMassVsPt->SetBinContent(i, sigma[i]*1000);
        sdtMassVsPt->SetBinError(i, sigmaErrors[i]*1000);
    }
    return sdtMassVsPt;
};



TH1D*  signalToNoise(ROOT::RDF::RResultPtr<TH2D> hMassVsPt, std::string resultName, std::string dir, std::string particle){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    std::unique_ptr<TFile> projections(TFile::Open(Form("projections/projectionsOmegaMass_%s.root",dir.c_str()), "RECREATE"));
 
    double pt[50];
    double mean[50];
    double signal[50];
    double noise[50];
    double ratioSN[50];
    double ratioSNError[50];

    float rangea, rangeb;

    if (particle == "Xi"){
        rangea = 1.300;
        rangeb = 1.345;
    }
    else if (particle == "Omega"){
        rangea = 1.637;
        rangeb = 1.707;
    }

    for (int i = 0; i < 50; i++){   

        TH1D *h = hMassVsPt->ProjectionY(Form("bin%d",i), i,i);
        h->GetXaxis()->SetRangeUser(rangea, rangeb);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.5);
    
        std::tuple<double,double,double,double> signal_backgroung_serr_berr = fitRoofit(h, rangea, rangeb, Form("StoN/bin%d.pdf", i), dir , "cb", particle, 0.1+ i*0.2, true, false);

        double signal = std::get<0>(signal_backgroung_serr_berr);
        double background = std::get<1>(signal_backgroung_serr_berr);
        double signalError = std::get<2>(signal_backgroung_serr_berr);
        double backgroundError = std::get<3>(signal_backgroung_serr_berr);


        std::cout<<"signal: "<<signal<<"  background: "<<background<<std::endl;

        pt[i] = 0.1 + (i)*0.2;     //////// !!!!!!!!!!!!!!!!!!!!
        ratioSN[i] = signal/(signal+background);
        ratioSNError[i]= sqrt(  pow(  (background/pow(signal+background,2))  * signalError     ,2) + pow( (-(signal/pow(signal+background,2)))  * signalError ,2)    );
        h->Write();
        h->SetStats(0);
        
    }
    projections->Close();


    TH1D * SToN = new TH1D("SToN", "", 50, 0, 10);

    TCanvas c;

    for (int i=6; i<36; i++){
        SToN->SetBinContent(i, ratioSN[i]);
        SToN->SetBinError(i,ratioSNError[i]);
    }
    std::string title = ("; #it{p}_{T} (GeV/c); #frac{S}{S+B}");
    SToN->SetTitle(title.c_str());
    // SToN->SetMarkerColor(8);
    SToN->SetMarkerStyle(21);
    SToN->SetName(resultName.c_str());
    SToN->GetXaxis()->SetRangeUser(0.9,7.1);
    SToN->GetYaxis()->SetRangeUser(0.5,1.5);
    return SToN;
};

TH1F* sigmaDCAvsPt(ROOT::RDF::RResultPtr<TH2D> hDCAvsPt, std::string resultName, std::string axis, std::string dir, std::string particle, std::string type){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // std::cout<<"------------"<<std::endl;

    std::unique_ptr<TFile> projections(TFile::Open(Form("projections/projections_%s.root",dir.c_str()), "RECREATE"));
    double sigma[25];
    double sigmaErrors[25];
    double mean[25];
    double meanErrors[25];
    double pt[25];


    int countbin = 0;
    for (int i = 0; i < 50; i=i+2){   /// skipping the first two bin (400 Mev)

        TH1D *h = hDCAvsPt->ProjectionY(Form("bin%d",i), i,i+1);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.5);

        std::tuple<double,double,double,double> mean_sigma_merr_serr = fitRoofit(h,-0.02, 0.02, Form("sigmaDCAvsPt/%sbin%d.pdf", type.c_str(),i), dir , "cb", particle, 0.2 + (countbin)*0.4,false, true, axis);

        // std::cout<<"\n+++++++++++++++++++++++++++sigma "<<std::get<1>(mean_sigma_merr_serr)<<" "<<std::get<3>(mean_sigma_merr_serr)<<" ++++++++++++++++"<<std::endl;

        sigma[countbin] = std::get<1>(mean_sigma_merr_serr);
        sigmaErrors[countbin] = std::get<3>(mean_sigma_merr_serr);
        mean[countbin] = std::get<0>(mean_sigma_merr_serr);
        meanErrors[countbin] = std::get<2>(mean_sigma_merr_serr);
        pt[countbin] = 0.2 + (countbin)*0.4;    ////////////////////// !!!!!!!!!!!!!!!!!
        // std::cout<<"i: "<<i<<" pt: "<<pt[countbin]<<std::endl;
        h->Write();
        

        h->SetStats(0);
        countbin++;
    }

    projections->Close();

    // TGraphErrors* DCAresolution = new TGraphErrors(24, pt, sigma, 0, sigmaErrors) ;
    TH1F* DCAresolution = new TH1F("DCAresolution", "", 25, 0, 10);

    // TCanvas c;
    // DCAresolution->Draw("AP");
    std::string title = Form("; #it{p}_{T} (GeV/c); #sigma_{DCA_{%s}} (#mum)", axis.c_str());
    DCAresolution->SetTitle(title.c_str());
    DCAresolution->SetMarkerStyle(21);
    DCAresolution->SetName(resultName.c_str());
    DCAresolution->GetXaxis()->SetRangeUser(0,10);
    DCAresolution->GetYaxis()->SetRangeUser(0,80);
  
    for (int i=0; i<25; i++){
        DCAresolution->SetBinContent(i+1, sigma[i]);
        DCAresolution->SetBinError(i+1, sigmaErrors[i]);
    }


    // std::unique_ptr<TGraphErrors> DCAresolution(50, pt, sigma, 0, sigmaErrors);
    // c.SaveAs(Form("%s.png", resultName.c_str()));

    // DCAresolution->Write();
    return DCAresolution;
};


// template<typename T>
void fillHistograms(ROOT::RDF::RNode table, TFile* myFile, std::string dir, std::string particle){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    /// histograms 1D: masses
    auto hMassV0 = table.Histo1D({"hMassV0", ";Inv. mass (GeV/#it{c}^{2});Counts", 125, 1.090, 1.140},"fMassV0");
    //v0
    fitRoofit(hMassV0.GetPtr(),1.090, 1.140, "MassV0.pdf", dir , "cb", "it{}V0", -1);
    
    if (particle == "Cascade"){
        auto hMassCascadevsPt1 = table.Histo2D({"hMassCascadevsPt1", ";#it{p}_{T} (GeV/#it{c});Inv. mass (GeV/#it{c}^{2})", 50, 0.,10., 1115, 1.310, 1.695}, "fCascPt", "fMassXi");
        auto hMassCascadevsPt2 = table.Histo2D({"hMassCascadevsPt2", ";#it{p}_{T} (GeV/#it{c});Inv. mass (GeV/#it{c}^{2})", 50, 0.,10., 1115, 1.310, 1.695}, "fCascPt", "fMassOmega");
        myFile->cd(dir.c_str());
        hMassCascadevsPt1->Write();
        hMassCascadevsPt2->Write();
    }

    else if (particle=="Xi" || particle == "Cascade"){

        auto hMassXi = table.Histo1D({"hMassXi", ";Inv. mass (GeV/#it{c}^{2});Counts", 125, 1.310, 1.335},"fMassXi");
        
        /// fitting masses:
        //xi
        std::tuple<double,double,double,double> mean_sigma_merr_serr = fitRoofit(hMassXi.GetPtr(),1.310, 1.335, "MassXi.pdf", dir , "cb", particle, -1);

        /// histograms 2D: mass vs DCA
        auto hXiMassVsDCAxy = table.Histo2D({"hXiMassVsDCAxy", ";#Xi DCA_{xy} (cm);Inv. mass (GeV/#it{c}^{2})", 100, -0.02, 0.02, 100, 1.300, 1.400}, "fCascDCAxy", "fMassXi");
        auto hXiMassVsDCAz = table.Histo2D({"hXiMassVsDCAz", ";#Xi DCA_{z} (cm);Inv. mass (GeV/#it{c}^{2})", 100, -0.02, 0.02, 100, 1.300, 1.400}, "fCascDCAz", "fMassXi");


        /// histograms 2D: DCA vs pt (mass cut: 2 sigma)
        // xi
        double meanMassXi = std::get<0>(mean_sigma_merr_serr);
        double sigmaMassXi = std::get<1>(mean_sigma_merr_serr);


     
        auto casctableCutXiMass = table.Filter([&](double fMassXi){if (fMassXi > meanMassXi - 2. * sigmaMassXi && fMassXi < meanMassXi + 2. * sigmaMassXi) return true; else return false;},{"fMassXi"});
        auto hXiCutsDCAxyVsPt = casctableCutXiMass.Histo2D({"hXiCutsDCAxyVsPt", ";#it{p}_{T} (GeV/#it{c});#Xi DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
        auto hXiCutsDCAzVsPt = casctableCutXiMass.Histo2D({"hXiCutsDCAzVsPt", ";#it{p}_{T} (GeV/#it{c});#Xi DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");
        
        auto hXiDCAxy = casctableCutXiMass.Histo1D({"hXiDCAxy", ";#Xi DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
        fitRoofit(hXiDCAxy.GetPtr(),-0.02, 0.02, "XiDCAxy.pdf", dir , "cb", particle, -1,false, true, "xy");
        auto hXiDCAz = casctableCutXiMass.Histo1D({"hXiDCAz", ";#Xi DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");        
        fitRoofit(hXiDCAz.GetPtr(),-0.02, 0.02, "XiDCAz.pdf", dir , "cb", particle, -1,false, true, "z");

        auto hXiDCAxyCuts = casctableCutXiMass.Histo1D({"hXiDCAxyCuts", ";#Xi DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
        fitRoofit(hXiDCAxyCuts.GetPtr(),-0.02, 0.02, "XiDCAxyCuts.pdf", dir , "cb", particle, -1,false, true, "xy");
        auto hXiDCAzCuts = casctableCutXiMass.Histo1D({"hXiDCAzCuts", ";#Xi DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");        
        fitRoofit(hXiDCAzCuts.GetPtr(),-0.02, 0.02, "XiDCAzCuts.pdf", dir , "cb", particle, -1,false, true, "z");

          
        TH1F * DCAxyResolutionXiCuts = sigmaDCAvsPt(hXiCutsDCAxyVsPt, "DCAxyResolutionXiCuts", "xy", dir.substr(0, dir.size()-1), particle, "cuts/xy_");
        TH1F * DCAzResolutionXiCuts = sigmaDCAvsPt(hXiCutsDCAzVsPt, "DCAzResolutionXiCuts","z", dir.substr(0, dir.size()-1), particle, "cuts/z_");

        auto casctableTailsXiMass = table.Filter([&](double fMassXi){if (fMassXi < meanMassXi - 4. * sigmaMassXi || fMassXi > meanMassXi + 4. * sigmaMassXi) return true; else return false;},{"fMassXi"});
        auto hXiTailsDCAxyVsPt = casctableTailsXiMass.Histo2D({"hXiTailsDCAxyVsPt", ";#it{p}_{T} (GeV/#it{c});#Xi DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
        auto hXiTailsDCAzVsPt = casctableTailsXiMass.Histo2D({"hXiTailsDCAzVsPt", ";#it{p}_{T} (GeV/#it{c});#Xi DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");

        auto hXiDCAxyTails = casctableTailsXiMass.Histo1D({"hXiDCAxyTails", ";#Xi DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
        fitRoofit(hXiDCAxyTails.GetPtr(),-0.02, 0.02, "XiDCAxyTails.pdf", dir , "cb", particle, -1,false, true, "xy");
        auto hXiDCAzTails = casctableTailsXiMass.Histo1D({"hXiDCAzTails", ";#Xi DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");        
        fitRoofit(hXiDCAzTails.GetPtr(),-0.02, 0.02, "XiDCAzTails.pdf", dir , "cb", particle, -1,false, true, "z");
       
        TH1F * DCAxyResolutionXiTails = sigmaDCAvsPt(hXiTailsDCAxyVsPt, "DCAxyResolutionXiTails", "xy", dir.substr(0, dir.size()-1), particle, "tails/xy_");
        TH1F * DCAzResolutionXiTails = sigmaDCAvsPt(hXiTailsDCAzVsPt, "DCAzResolutionXiTails","z", dir.substr(0, dir.size()-1), particle,"tails/z_");

        // histograms 2D: mass vs pt
        auto hXiMassVsPt = table.Histo2D({"hXiMassVsPt", ";#it{p}_{T} (GeV/#it{c});Inv. mass (GeV/#it{c}^{2})", 50, 0.,10., 125, 1.300, 1.345}, "fCascPt", "fMassXi");

        //extract sigma(mass) vs pt for each mass
        TH1F * stdMassXi = stdMass(hXiMassVsPt, "stdMassXi", dir.substr(0, dir.size()-1), "Xi", "; #it{p}_{T} (GeV/#it{c}); #sigma_{Mass_{#Xi}} (MeV/#it{c}^{2})");

         //signal/noise 
        TH1D* StoN = signalToNoise(hXiMassVsPt, "StoN", dir.substr(0, dir.size()-1), "Xi");

        myFile->cd(dir.c_str());


        StoN->Write();
        hMassXi->Write();
        hXiMassVsDCAxy->Write();
        hXiMassVsDCAz->Write();
        hXiCutsDCAxyVsPt->Write();
        hXiCutsDCAzVsPt->Write();
        hXiDCAxy->Write();
        hXiDCAz->Write();
        hXiDCAxyCuts->Write();
        hXiDCAzCuts->Write();
        hXiDCAxyTails->Write();
        hXiDCAzTails->Write();
        DCAxyResolutionXiCuts->Write();
        DCAzResolutionXiCuts->Write();
        hXiTailsDCAxyVsPt->Write();
        hXiTailsDCAzVsPt->Write();
        DCAxyResolutionXiTails->Write();
        DCAzResolutionXiTails->Write();
        hXiMassVsPt->Write();
        stdMassXi->Write();


        TCanvas c0;
        // StoN->GetYaxis()->SetRangeUser(-0.5, 0.5);
        TF1 fitStoN = TF1("fitStoN", "pol0", 0.5, 7.5);
        fitStoN.SetLineWidth(1);
        StoN->Fit("fitStoN", "R");
        auto p = fitStoN.GetParameter(0);
        auto perr = fitStoN.GetParError(0);
        StoN->SetMarkerSize(0.5);
        StoN->SetMarkerStyle(20);
        StoN->Draw();
        fitStoN.Draw("same");
        TPaveText *pv = new TPaveText(0.6, 0.7, 0.8, 0.9, "NDC");
        pv->SetTextFont(132);
        pv->SetFillStyle(0);
        pv->SetBorderSize(0);
        pv->AddText(Form("#%s", particle.c_str()));
        ((TText*)pv->GetListOfLines()->Last())->SetTextSize(0.07);
       
        pv->AddText(Form("Avg: %.3f #pm %.3f",  p, perr));
        pv->Draw("same");

        c0.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/StoN.pdf").c_str());



        TCanvas c;
        // stdMassXi->GetYaxis()->SetRangeUser(0., 0.003);
        stdMassXi->SetMarkerSize(0.5);
        stdMassXi->SetMarkerStyle(20);
        stdMassXi->Draw();
        c.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/stdMassXi.pdf").c_str());

        TCanvas c1;
        DCAxyResolutionXiCuts->SetMarkerSize(0.5);
        DCAxyResolutionXiCuts->SetMarkerStyle(20);
        DCAxyResolutionXiCuts->Draw();
        c1.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAxyResolutionXiCuts.pdf").c_str());

        TCanvas c2;
        DCAzResolutionXiCuts->SetMarkerSize(0.5);
        DCAzResolutionXiCuts->SetMarkerStyle(20);
        DCAzResolutionXiCuts->Draw();
        c2.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAzResolutionXiCuts.pdf").c_str());


        TCanvas c3;
        DCAxyResolutionXiTails->SetMarkerSize(0.5);
        DCAxyResolutionXiTails->SetMarkerStyle(20);
        DCAxyResolutionXiTails->Draw();
        c3.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAxyResolutionXiTails.pdf").c_str());

        TCanvas c4;
        DCAzResolutionXiTails->SetMarkerSize(0.5);
        DCAzResolutionXiTails->SetMarkerStyle(20);
        DCAzResolutionXiTails->Draw();
        c4.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAzResolutionXiTails.pdf").c_str());

        std::cout<<"++++++++++++++++++++++++++meanMassXi: "<<meanMassXi<<" sigmaMassXi: "<<sigmaMassXi<<std::endl;



    }


    else if (particle == "Omega" || particle == "Cascade"){

        auto hMassOmega = table.Histo1D({"hMassOmega", ";Inv. mass (GeV/#it{c}^{2});Counts", 125, 1.600, 1.750},"fMassOmega");
        std::tuple<double, double,double,double> mean_sigma_merr_serr = fitRoofit(hMassOmega.GetPtr(),1.647, 1.697, "MassOmega.pdf", dir , "cb", particle, -1);

        /// histograms 2D: mass vs DCA
        auto hOmegaMassVsDCAxy = table.Histo2D({"hOmegaMassVsDCAxy", ";#Omega DCA_{xy} (cm);Inv. mass (GeV/#it{c}^{2})", 100, -0.02, 0.02, 100, 1.600, 1.800}, "fCascDCAxy", "fMassOmega");
        auto hOmegaMassVsDCAz = table.Histo2D({"hOmegaMassVsDCAz", ";#Omega DCA_{z} (cm);Inv. mass (GeV/#it{c}^{2})", 100, -0.02, 0.02, 100, 1.600, 1.800}, "fCascDCAz", "fMassOmega");

        /// histograms 2D: DCA vs pt (mass cut: 2 sigma)
        //omega
        double meanMassOmega = std::get<0>(mean_sigma_merr_serr);
        double sigmaMassOmega = std::get<1>(mean_sigma_merr_serr);

        auto casctableCutOmegaMass = table.Filter([&](double fMassOmega){if (fMassOmega > meanMassOmega - 2. * sigmaMassOmega && fMassOmega < meanMassOmega + 2. * sigmaMassOmega) return true; else return false;},{"fMassOmega"});
        auto hOmegaCutsDCAxyVsPt = casctableCutOmegaMass.Histo2D({"hOmegaCutsDCAxyVsPt", ";#it{p}_{T} (GeV/#it{c});#Omega DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
        auto hOmegaCutsDCAzVsPt = casctableCutOmegaMass.Histo2D({"hOmegaCutsDCAzVsPt", ";#it{p}_{T} (GeV/#it{c});#Omega DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");

        auto hOmegaDCAxy = casctableCutOmegaMass.Histo1D({"hOmegaDCAxy", ";#Omega DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
        fitRoofit(hOmegaDCAxy.GetPtr(),-0.02, 0.02, "OmegaDCAxy.pdf", dir , "cb", particle, -1, false, true, "xy");
        auto hOmegaDCAz = casctableCutOmegaMass.Histo1D({"hOmegaDCAz", ";#Omega DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");
        fitRoofit(hOmegaDCAz.GetPtr(),-0.02, 0.02, "OmegaDCAz.pdf", dir , "cb", particle, -1,false, true, "z");
        
        auto hOmegaDCAxyCuts = casctableCutOmegaMass.Histo1D({"hOmegaDCAxyCuts", ";#Omega DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
        fitRoofit(hOmegaDCAxyCuts.GetPtr(),-0.02, 0.02, "OmegaDCAxyCuts.pdf", dir , "cb", particle, -1,false, true, "xy");
        auto hOmegaDCAzCuts = casctableCutOmegaMass.Histo1D({"hOmegaDCAzCuts", ";#Omega DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");        
        fitRoofit(hOmegaDCAzCuts.GetPtr(),-0.02, 0.02, "OmegaDCAzCuts.pdf", dir , "cb", particle, -1,false, true, "z");



        TH1F * DCAxyResolutionOmegaCuts = sigmaDCAvsPt(hOmegaCutsDCAxyVsPt, "DCAxyResolutionOmegaCuts", "xy", dir.substr(0, dir.size()-1), particle, "cuts/xy_");
        TH1F * DCAzResolutionOmegaCuts = sigmaDCAvsPt(hOmegaCutsDCAzVsPt, "DCAzResolutionOmegaCuts","z", dir.substr(0, dir.size()-1), particle, "cuts/z_");

        auto casctableTailsOmegaMass = table.Filter([&](double fMassOmega){if (fMassOmega < meanMassOmega -4. * sigmaMassOmega || fMassOmega > meanMassOmega + 4. * sigmaMassOmega) return true; else return false;},{"fMassOmega"});
        auto hOmegaTailsDCAxyVsPt = casctableTailsOmegaMass.Histo2D({"hOmegaTailsDCAxyVsPt", ";#it{p}_{T} (GeV/#it{c});#Omega DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
        auto hOmegaTailsDCAzVsPt = casctableTailsOmegaMass.Histo2D({"hOmegaTailsDCAzVsPt", ";#it{p}_{T} (GeV/#it{c});#Omega DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");

        auto hOmegaDCAxyTails = casctableTailsOmegaMass.Histo1D({"hOmegaDCAxyTails", ";#Omega DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
        fitRoofit(hOmegaDCAxyTails.GetPtr(),-0.02, 0.02, "OmegaDCAxyTails.pdf", dir , "cb", particle, -1,false, true, "xy");
        auto hOmegaDCAzTails = casctableTailsOmegaMass.Histo1D({"hOmegaDCAzTails", ";#Omega DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");        
        fitRoofit(hOmegaDCAzTails.GetPtr(),-0.02, 0.02, "OmegaDCAzTails.pdf", dir , "cb", particle, -1,false, true, "z");
     
        TH1F * DCAxyResolutionOmegaTails = sigmaDCAvsPt(hOmegaTailsDCAxyVsPt, "DCAxyResolutionOmegaTails", "xy", dir.substr(0, dir.size()-1), particle, "tails/xy_");
        TH1F * DCAzResolutionOmegaTails = sigmaDCAvsPt(hOmegaTailsDCAzVsPt, "DCAzResolutionOmegaTails","z", dir.substr(0, dir.size()-1), particle, "tails/z_");

        // histograms 2D: mass vs pt
        auto hOmegaMassVsPt = table.Histo2D({"hOmegaMassVsPt", ";#it{p}_{T} (GeV/#it{c});Inv. mass (GeV/#it{c}^{2})", 50, 0.,10., 100, 1.637, 1.707}, "fCascPt", "fMassOmega");
        // auto hOmegaMassVsPt2sigma = casctableCutOmegaMass.Histo2D({"hOmegaMassVsPt2sigma", ";#it{p}_{T} (GeV/#it{c});Inv. mass (GeV/#it{c}^{2})", 50, 0.,10., 100, 1.647, 1.697}, "fCascPt", "fMassOmega");
        
        // auto hOmegaMass2sigma = casctableCutOmegaMass.Histo1D({"hOmegaMass2sigma", " ;Inv. mass (GeV/#it{c}^{2}); ", 100, 1.647, 1.697}, "fMassOmega");

        //signal/noise 
        TH1D* StoN = signalToNoise(hOmegaMassVsPt, "StoN", dir.substr(0, dir.size()-1), "Omega");

        //extract sigma(mass) vs pt for each mass
        TH1F * stdMassOmega = stdMass(hOmegaMassVsPt, "stdMassOmega", dir.substr(0, dir.size()-1), "Omega", "; #it{p}_{T} (GeV/c); #sigma_{Mass_{#Omega}} (MeV/#it{c}^{2})");

        myFile->cd(dir.c_str());

        // hOmegaMass2sigma->Write();
        hMassOmega->Write();
        hOmegaMassVsDCAxy->Write();
        hOmegaMassVsDCAz->Write();
        hOmegaCutsDCAxyVsPt->Write();
        hOmegaCutsDCAzVsPt->Write();
        hOmegaDCAxy->Write();
        hOmegaDCAz->Write();
        DCAxyResolutionOmegaCuts->Write();
        DCAzResolutionOmegaCuts->Write();
        hOmegaTailsDCAxyVsPt->Write();
        hOmegaTailsDCAzVsPt->Write();
        DCAxyResolutionOmegaTails->Write();
        DCAzResolutionOmegaTails->Write();
        hOmegaMassVsPt->Write();
        StoN->Write();
        stdMassOmega->Write();

        TCanvas c;
        // stdMassOmega->GetYaxis()->SetRangeUser(0., 0.003);
        stdMassOmega->SetMarkerSize(0.5);
        stdMassOmega->SetMarkerStyle(20);
        stdMassOmega->Draw();

        std::cout<<"--------dir "<<dir<<"---------"<<std::endl;
        c.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/stdMassOmega.pdf").c_str());


        TCanvas c1;
        // StoN->GetYaxis()->SetRangeUser(-0.5, 0.5);
        TF1 fitStoN = TF1("fitStoN", "pol0", 0.5, 7.5);
        fitStoN.SetLineWidth(1);
        StoN->Fit("fitStoN");
        auto p = fitStoN.GetParameter(0);
        auto perr = fitStoN.GetParError(0);
        StoN->SetMarkerSize(0.5);
        StoN->SetMarkerStyle(20);
        StoN->Draw();
        fitStoN.Draw("same");
        TPaveText *pv = new TPaveText(0.6, 0.7, 0.8, 0.9, "NDC");
        pv->SetTextFont(132);
        pv->SetFillStyle(0);
        pv->SetBorderSize(0);
        pv->AddText(Form("#%s", particle.c_str()));
        ((TText*)pv->GetListOfLines()->Last())->SetTextSize(0.07);
       
        pv->AddText(Form("Avg: %.3f #pm %.3f",  p, perr));
        pv->Draw("same");

        c1.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/StoN.pdf").c_str());

        TCanvas c2; 
        DCAxyResolutionOmegaCuts->SetMarkerSize(0.5);
        DCAxyResolutionOmegaCuts->SetMarkerStyle(20);
        DCAxyResolutionOmegaCuts->Draw();
        c2.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAxyResolutionOmegaCuts.pdf").c_str());

        TCanvas c3;
        DCAzResolutionOmegaCuts->SetMarkerSize(0.5);
        DCAzResolutionOmegaCuts->SetMarkerStyle(20);
        DCAzResolutionOmegaCuts->Draw();
        c3.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAzResolutionOmegaCuts.pdf").c_str());

        TCanvas c4;
        DCAxyResolutionOmegaTails->SetMarkerSize(0.5);
        DCAxyResolutionOmegaTails->SetMarkerStyle(20);
        DCAxyResolutionOmegaTails->Draw();
        c4.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAxyResolutionOmegaTails.pdf").c_str());

        TCanvas c5;
        DCAzResolutionOmegaTails->SetMarkerSize(0.5);
        DCAzResolutionOmegaTails->SetMarkerStyle(20);
        DCAzResolutionOmegaTails->Draw();
        c5.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAzResolutionOmegaTails.pdf").c_str());



    }
    
    

    /// histograms 1D: DCA
    auto hCascDCAxy = table.Histo1D({"hCascDCAxy", ";Cascade DCA_{xy} (cm);Counts", 200, -0.02, 0.02},"fCascDCAxy");
    auto hCascDCAz = table.Histo1D({"hCascDCAz", ";Cascade DCA_{z} (cm);Counts", 200, -0.02, 0.02},"fCascDCAz");
    auto hProtonDCAxy = table.Histo1D({"hProtonDCAxy", ";Proton DCA_{xy} (cm);Counts", 200, -5., 5.},"fProtonDCAxy");
    auto hProtonDCAz = table.Histo1D({"hProtonDCAz", ";Proton DCA_{z} (cm);Counts", 200, -5., 5.},"fProtonDCAz");
    auto hPionDCAxy = table.Histo1D({"hPionDCAxy", ";Pion DCA_{xy} (cm);Counts", 200, -10., 10.},"fPionDCAxy");
    auto hPionDCAz = table.Histo1D({"hPionDCAz", ";Pion DCA_{z} (cm);Counts", 200, -10., 10.},"fPionDCAz");
    auto hBachDCAxy = table.Histo1D({"hBachDCAxy", ";Bachelor DCA_{xy} (cm);Counts", 200, -10., 10.},"fBachDCAxy");
    auto hBachDCAz = table.Histo1D({"hBachDCAz", ";Bachelor DCA_{z} (cm);Counts", 200, -10., 10.},"fBachDCAz");


    /// histograms 2D: DCA vs pt (no cuts)
    auto hDCAxyVsPt = table.Histo2D({"hDCAxyVsPt", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAxy");
    auto hDCAzVsPt = table.Histo2D({"hDCAzVsPt", ";#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", 50, 0.,10., 100, -0.02, 0.02}, "fCascPt", "fCascDCAz");
  
    TH1F * DCAxyResolution = sigmaDCAvsPt(hDCAxyVsPt, "DCAxyResolution", "xy", dir.substr(0, dir.size()-1),particle,"xy_");
    TH1F * DCAzResolution = sigmaDCAvsPt(hDCAzVsPt, "DCAzResolution","z", dir.substr(0, dir.size()-1), particle,"z_");
 
    myFile->cd(dir.c_str());

        
    hMassV0->Write();

    hCascDCAxy->Write();
    hCascDCAz->Write();
    hProtonDCAxy->Write();
    hProtonDCAz->Write();
    hPionDCAxy->Write();
    hPionDCAz->Write();
    hBachDCAxy->Write();
    hBachDCAz->Write();

    


    hDCAxyVsPt->Write();
    hDCAzVsPt->Write();

    DCAxyResolution->Write();
    DCAzResolution->Write();

    TCanvas c6;
    DCAxyResolution->SetMarkerSize(0.5);
    DCAxyResolution->SetMarkerStyle(20);
    DCAxyResolution->Draw();
    c6.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAxyResolution.pdf").c_str());

    TCanvas c7;
    DCAzResolution->SetMarkerSize(0.5);
    DCAzResolution->SetMarkerStyle(20);
    DCAzResolution->Draw();
    c7.SaveAs(("plots/" + dir.substr(0, dir.size()-1) + "/DCAzResolution.pdf").c_str());



}





void readAO2D(){

    gErrorIgnoreLevel = kWarning;
    gStyle->SetImageScaling(1.5);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);

    // Margins
    gStyle->SetPadBottomMargin(0.1); 
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.15);
    //Frame options
    gStyle->SetTitleSize(0.04,"xyz"); // size of axis title font
    gStyle->SetTitleSize(0.05,"");
    // gStyle->SetTitleFont(132,""); // font option
    // gStyle->SetTitleFont(132,"xyz"); // font option
    gStyle->SetTitleOffset(1.2, "yx");
    //Label options
    gStyle->SetLabelSize(0.035,"xyz"); // size of axis value font
    // gStyle->SetLabelFont(132,"xyz");
    //Legend options
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(132);
    gStyle->SetLegendFillColor(-1);
    //Stat options
    // gStyle->SetStatFont(132);
    // gStyle->SetStatFontSize(0.05);
    //Other
    gStyle->SetFillStyle(4000);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);


    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    int minClusTPC = 70;
    float minNsigmaTPC = -4.;
    float maxNsigmaTPC = 4.;

    ROOT::RDataFrame casctable("DF_2316989016393920/O2npcasctable", "AO2D_apass2.root");

    auto omegatable = casctable.Filter("abs(fMassXi-1.32171) > 0.005"); //(massXi - 1.32171) > 0.005 
    auto xitable = casctable.Filter("abs(fMassOmega-1.67245) > 0.005"); //(massXi - 1.67245) > 0.005
    auto casctablePIDcuts = casctable.Filter("fProtonNClusTPC>70 and fPionNClusTPC>70 and fProtonTPCNSigma > -4 and fProtonTPCNSigma<4 and fPionTPCNSigma>-4 and fPionTPCNSigma<4 and fBachKaonTPCNSigma>-4 and fBachKaonTPCNSigma<4 and fBachPionTPCNSigma>-4 and fBachPionTPCNSigma<4");
    auto omegatablePIDcuts = omegatable.Filter("fProtonNClusTPC>70 and fPionNClusTPC>70 and fProtonTPCNSigma > -4 and fProtonTPCNSigma<4 and fPionTPCNSigma>-4 and fPionTPCNSigma<4 and fBachKaonTPCNSigma>-4 and fBachKaonTPCNSigma<4");
    auto xitablePIDcuts = casctable.Filter("fProtonNClusTPC>70 and fPionNClusTPC>70 and fProtonTPCNSigma > -4 and fProtonTPCNSigma<4 and fPionTPCNSigma>-4 and fPionTPCNSigma<4 and fBachKaonTPCNSigma>-4 and fBachKaonTPCNSigma<4 and fBachPionTPCNSigma>-4 and fBachPionTPCNSigma<4");
  
    TFile* myFile(TFile::Open("histos.root", "RECREATE"));
    myFile->mkdir("Cascade/");
    myFile->mkdir("Omega/");
    myFile->mkdir("Xi/");
    myFile->mkdir("Cascade+PIDcuts/");
    myFile->mkdir("Omega+PIDcuts/");
    myFile->mkdir("Xi+PIDcuts/");
    
    // fillHistograms(casctable, myFile, "Cascade/", "Cascade");
    // fillHistograms(omegatable, myFile, "Omega/", "Omega");
    // fillHistograms(xitable, myFile, "Xi/", "Xi");
    // fillHistograms(casctablePIDcuts, myFile, "Cascade+PIDcuts/", "Cascade");
    fillHistograms(omegatablePIDcuts, myFile, "Omega+PIDcuts/", "Omega");
    fillHistograms(xitablePIDcuts, myFile, "Xi+PIDcuts/", "Xi");

    myFile->Close();
}