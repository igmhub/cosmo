// Created 7-Feb-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macro to plot results of the baofit program.

#include <fstream>
#include <iostream>

void plotBaoFit(const char *filename = "fit.dat") {
    // Initialize graphics options.
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    Double_t rGradient[3] = { 0.0, 1.0, 1.0 };
    Double_t gGradient[3] = { 0.0, 1.0, 0.0 };
    Double_t bGradient[3] = { 1.0, 1.0, 0.0 };
    Double_t sGradient[3] = { 0.0, 0.5, 1.0 };
    TColor::CreateGradientColorTable(3,sGradient,rGradient,gGradient,bGradient,256);
    //gStyle->SetPalette(1,0);
    
    // Sound horizon in Mpc/h at zdrag = 1020.49, calculated from Eisenstein & Hu 1997 using the command:
    // cosmocalc --omega-lambda 0.734 --omega-matter 0.266 --omega-baryon 0.0449 --hubble-constant 0.710
    Double_t r3dContours[1] = { 108.719 };
    
    std::ifstream in(filename);
    // Read the binning parameters.
    int nll,nsep,nz,ndata,oversampling,ncontour;
    double minll,minsep,minz,dll,dsep,dz,scale;
    in >> nll >> minll >> dll;
    in >> nsep >> minsep >> dsep;
    in >> nz >> minz >> dz;
    in >> ndata >> oversampling >> ncontour >> scale;
    
    r3dContours[0] *= scale;
    
    // Initialize the drawing canvas.
    TCanvas *canvas = new TCanvas("canvas","canvas",nz*400,800);
    canvas->SetMargin(0.05,0.01,0.05,0.05);
    canvas->Divide(nz,2,0.,0.);

    // Create 2D histograms in (ll,sep) for each redshift.
    double cLight = 299.792458; // 1e3 km/s
    double vmin = cLight*minll, vmax = cLight*(minll+nll*dll);
    for(int iz = 0; iz < nz; ++iz) {
        double z(minz + (iz+0.5)*dz);
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        if(dataHist) continue;
        TH2F *dataHist = new TH2F(Form("data%d",iz),Form("Fit data for z = %.2f",z),
            nsep,minsep,minsep+nsep*dsep,nll,vmin,vmax);
        dataHist->SetXTitle("Angular separation (arcmin)");
        dataHist->SetYTitle("Relative radial velocity (10^{3}km/s)");
        TH2F *pullHist = dataHist->Clone(Form("pull%d",iz));
        TH2F *modelHist = new TH2F(Form("model%d",iz),Form("Model predictions for z = %.2f",z),
            oversampling*nsep,minsep,minsep+nsep*dsep,oversampling*nll,vmin,vmax);
        TH2F *r3dHist = new TH2F(Form("r3d%d",iz),Form("Co-moving 3D separation for z = %.2f",z),
            oversampling*nsep,minsep,minsep+nsep*dsep,oversampling*nll,vmin,vmax);
    }
    
    // Loop over bin data in the input file.
    int index;
    double data,pull;
    double *sumSqData = new double[nz];
    for(int iz = 0; iz < nz; ++iz) {
        sumSqData[iz] = 0;
    }
    for(int k = 0; k < ndata; ++k) {
        in >> index >> data >> pull;
        int iz = index % nz;
        sumSqData[iz] += data*data;
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        TH2F *pullHist = (TH2F *)gDirectory->Get(Form("pull%d",iz));
        int isep = 1 + (index/nz) % nsep;
        int ill = 1 + (index/(nz*nsep)) % nll;
        dataHist->SetBinContent(isep,ill,data);
        pullHist->SetBinContent(isep,ill,pull);
    }
    for(int iz = 0; iz < nz; ++iz) {
        std::cout << "RMS = " << std::sqrt(sumSqData[iz]) << std::endl;
    }
    
    // Loop over model predictions in the input file.
    double r3d,pred;
    for(int iz = 0; iz < nz; ++iz) {
        TH2F *modelHist = (TH2F *)gDirectory->Get(Form("model%d",iz));
        TH2F *r3dHist = (TH2F *)gDirectory->Get(Form("r3d%d",iz));
        for(int isep = 0; isep < oversampling*nsep; ++isep) {
            for(int ill = 0; ill < oversampling*nll; ++ill) {
                in >> r3d >> pred;
                r3dHist->SetBinContent(isep+1,ill+1,r3d);
                modelHist->SetBinContent(isep+1,ill+1,pred);
            }
        }
    }

    // Draw plots.
    double nsig=3;
    for(iz = 0; iz < nz; ++iz) {
        double zmax = 0.05*std::sqrt(sumSqData[iz]);
        TH2F *dataHist = (TH2F *)gDirectory->Get(Form("data%d",iz));
        TH2F *pullHist = (TH2F *)gDirectory->Get(Form("pull%d",iz));
        TH2F *modelHist = (TH2F *)gDirectory->Get(Form("model%d",iz));
        TH2F *r3dHist = (TH2F *)gDirectory->Get(Form("r3d%d",iz));

        canvas->cd(iz+1);
        canvas->GetPad(iz+1)->SetMargin(0.11,0.03,0.10,0.01);
        dataHist->GetYaxis()->SetTitleOffset(1.3);
        dataHist->SetMaximum(+zmax);
        dataHist->SetMinimum(-zmax*0.99); // factor of 0.99 ensures that zero is white
        // Trunctate bins outside the limits so that they are colored correctly.
        for(int ill = 0; ill < nll; ++ill) {
            for(int isep = 0; isep < nsep; ++isep) {
                double data = dataHist->GetBinContent(isep+1,ill+1);
                if(data < -zmax) data = -0.99*zmax;
                if(data > +zmax) data = +0.99*zmax;
                dataHist->SetBinContent(isep+1,ill+1,data);
                double pull = pullHist->GetBinContent(isep+1,ill+1);
                if(pull < -nsig) pull = -0.98*nsig;
                if(pull > +nsig) pull = +0.98*nsig;
                pullHist->SetBinContent(isep+1,ill+1,pull);
            }
        }
        dataHist->Draw("col");
        modelHist->SetMaximum(+zmax);
        modelHist->SetMinimum(-zmax);
        modelHist->Draw("cont3same");
        r3dHist->SetContour(1,r3dContours);
        r3dHist->SetLineWidth(5);
        r3dHist->SetLineColor(kGreen-6);
        r3dHist->Draw("cont3same");

        canvas->cd(nz+iz+1);
        canvas->GetPad(nz+iz+1)->SetMargin(0.11,0.03,0.10,0.01);
        pullHist->GetYaxis()->SetTitleOffset(1.3);
        pullHist->SetMaximum(+nsig);
        pullHist->SetMinimum(-nsig*0.99); // factor of 0.99 ensures that zero is white
        pullHist->Draw("col");
    }

    if(ncontour > 0) {
        TCanvas *canvas2 = new TCanvas("canvas2","canvas2",800,800);
        canvas2->Divide(2,2,0.001,0.001);
        Double_t *xContour = new Double_t[ncontour+1], *yContour = new Double_t[ncontour+1];
        TGraph *contourGraph[8];
        for(int ig = 0; ig < 8; ++ig) {
            for(int i = 0; i < ncontour; ++i) in >> xContour[i] >> yContour[i];
            // Close the contour
            xContour[ncontour] = xContour[0];
            yContour[ncontour] = yContour[0];
            contourGraph[ig] = new TGraph(ncontour+1,xContour,yContour);
            int ipad = ig%4+1;
            canvas2->cd(ipad);
            canvas2->GetPad(ipad)->SetMargin(0.15,0.02,0.10,0.02);
            canvas2->GetPad(ipad)->SetGridx();
            canvas2->GetPad(ipad)->SetGridy();
            contourGraph[ig]->SetLineColor(kBlue-8);
            contourGraph[ig]->SetLineWidth(3);
            contourGraph[ig]->Draw(ig < 4 ? "ALW":"L");
        }
        contourGraph[0]->GetHistogram()->SetXTitle("BAO Relative Scale");
        contourGraph[0]->GetHistogram()->SetYTitle("BAO Relative Amplitude");
        contourGraph[0]->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
        contourGraph[1]->GetHistogram()->SetXTitle("Lyman-alpha Bias");
        contourGraph[1]->GetHistogram()->SetYTitle("BAO Relative Amplitude");
        contourGraph[1]->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
        contourGraph[2]->GetHistogram()->SetXTitle("BAO Relative Scale");
        contourGraph[2]->GetHistogram()->SetYTitle("Redshift Distortion Beta");
        contourGraph[2]->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
        contourGraph[3]->GetHistogram()->SetXTitle("Lyman-alpha Bias");
        contourGraph[3]->GetHistogram()->SetYTitle("Redshift Distortion Beta");
        contourGraph[3]->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    }

    in.close();
}
