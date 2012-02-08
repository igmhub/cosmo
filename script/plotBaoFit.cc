// Created 7-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macro to plot results of the baofit program.

#include <fstream>
#include <iostream>

void plotBaoFit(const char *filename = "fit.dat") {
    // Initialize graphics options.
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    Double_t r[3] = { 0.0, 1.0, 1.0 };
    Double_t g[3] = { 0.0, 1.0, 0.0 };
    Double_t b[3] = { 1.0, 1.0, 0.0 };
    Double_t s[3] = { 0.0, 0.5, 1.0 };
    TColor::CreateGradientColorTable(3,s,r,g,b,256);
    //gStyle->SetPalette(1,0);
    
    std::ifstream in(filename);
    // Read the binning parameters.
    int nll,nsep,nz,ndata,oversampling;
    double minll,minsep,minz,dll,dsep,dz;
    in >> nll >> minll >> dll;
    in >> nsep >> minsep >> dsep;
    in >> nz >> minz >> dz;
    in >> ndata >> oversampling;
    
    // Initialize the drawing canvas.
    canvas = new TCanvas("canvas","canvas",nz*400,800);
    canvas->Divide(nz,2);

    // Create 2D histograms in (ll,sep) for each redshift.
    for(int iz = 0; iz < nz; ++iz) {
        double z(minz + (iz+0.5)*dz);
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        if(dataHist) continue;
        TH2F *dataHist = new TH2F(Form("data%d",iz),Form("Fit data for z = %.2f",z),
            nsep,minsep,minsep+nsep*dsep,nll,minll,minll+nll*dll);
        dataHist->SetXTitle("Pair angular separation (arcmin)");
        dataHist->SetYTitle("Pair LOS separation log(lam2/lam1)");
        TH2F *pullHist = dataHist->Clone(Form("pull%d",iz));
        TH2F *modelHist = new TH2F(Form("model%d",iz),Form("Model predictions for z = %.2f",z),
            oversampling*nsep,minsep,minsep+nsep*dsep,oversampling*nll,minll,minll+nll*dll);
    }
    
    // Loop over bin data in the input file.
    int index;
    double data,pull;
    for(int k = 0; k < ndata; ++k) {
        in >> index >> data >> pull;
        int iz = index % nz;
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        TH2F *pullHist = (TH2F *)gDirectory->Get(Form("pull%d",iz));
        int isep = 1 + (index/nz) % nsep;
        int ill = 1 + (index/(nz*nsep)) % nll;
        dataHist->SetBinContent(isep,ill,data);
        pullHist->SetBinContent(isep,ill,pull);
    }
    
    // Loop over model predictions in the input file.
    double pred;
    for(int iz = 0; iz < nz; ++iz) {
        TH2F *modelHist = (TH2F *)gDirectory->Get(Form("model%d",iz));
        for(int isep = 0; isep < oversampling*nsep; ++isep) {
            for(int ill = 0; ill < oversampling*nll; ++ill) {
                in >> pred;
                modelHist->SetBinContent(isep+1,ill+1,pred);
            }
        }
    }

    // Draw plots.
    double zmax=1e-4;
    for(iz = 0; iz < nz; ++iz) {
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        TH2F *modelHist = (TH2F *)gDirectory->Get(Form("model%d",iz));
        canvas->cd(iz+1);
        //canvas->GetPad(iz+1)->SetLogz();
        dataHist->SetMaximum(+zmax);
        dataHist->SetMinimum(-zmax);
        dataHist->Draw("col");
        modelHist->SetMaximum(+zmax);
        modelHist->SetMinimum(-zmax);
        modelHist->Draw("cont3same");
        TH2F *pullHist = (TH2F *)gDirectory->Get(Form("pull%d",iz));
        pullHist->SetMinimum(-3);
        pullHist->SetMaximum(+3);
        canvas->cd(nz+iz+1);
        pullHist->Draw("col");                
    }    

    in.close();
}