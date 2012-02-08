// Created 7-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macro to plot results of the baofit program.

#include <fstream>
#include <iostream>

void plotBaoFit(const char *filename = "fit.dat") {
    // Initialize graphics options.
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(1,0);
    
    std::ifstream in(filename);
    // Read the binning parameters.
    int nll,nsep,nz;
    double minll,minsep,minz,dll,dsep,dz;
    in >> nll >> minll >> dll;
    in >> nsep >> minsep >> dsep;
    in >> nz >> minz >> dz;
    
    // Initialize the drawing canvas.
    canvas = new TCanvas("canvas","canvas",nz*400,800);
    canvas->Divide(nz,2);

    // Create 2D histograms in (ll,sep) for each redshift.
    for(iz = 0; iz < nz; ++iz) {
        double z(minz + (iz+0.5)*dz);
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        if(dataHist) continue;
        TH2F *dataHist = new TH2F(Form("data%d",iz),Form("Fit data for z = %.2f",z),
            nsep,minsep,minsep+nsep*dsep,nll,minll,minll+nll*dll);
        dataHist->SetXTitle("Pair angular separation (arcmin)");
        dataHist->SetYTitle("Pair LOS separation log(lam2/lam1)");
        TH2F *pullHist = dataHist->Clone(Form("pull%d",iz));
    }
    
    // Loop over bin data in the input file.
    while(in.good() && !in.eof()) {
        int index;
        double data,pull;
        in >> index >> data >> pull;
        std::cout << index << ' ' << data << ' ' << pull << std::endl;
        if(in.eof()) break;
        if(!in.good()) {
            std::cout << "Error reading input file." << std::endl;
            return;
        }
        int iz = index % nz;
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        TH2F *pullHist = (TH2F *)gDirectory->Get(Form("pull%d",iz));
        int isep = 1 + (index/nz) % nsep;
        int ill = 1 + (index/(nz*nsep)) % nll;
        dataHist->SetBinContent(isep,ill,data);
        pullHist->SetBinContent(isep,ill,pull);
    }

    // Draw plots.
    for(iz = 0; iz < nz; ++iz) {
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        canvas->cd(iz+1);
        canvas->GetPad(iz+1)->SetLogz();
        dataHist->Draw("col");
        TH2F *pullHist = (TH2F *)gDirectory->Get(Form("pull%d",iz));
        pullHist->SetMinimum(-2);
        pullHist->SetMaximum(+2);
        canvas->cd(nz+iz+1);
        pullHist->Draw("col");                
    }    

    in.close();
}