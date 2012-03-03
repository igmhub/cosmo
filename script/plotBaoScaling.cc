// Created 2-Mar-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macro to plot the scaling of BAO parameter errors with bootstrap sample size.

TCanvas *canvas;
std::vector<double> sampleSize,scaleError,ampError;

void analyze(int nbs, int color) {
    TFile *file = new TFile(Form("results/bs%d.root",nbs),"READ");
    TTree *tree = (TTree*)file->Get("bs");

    bool first = (0 == sampleSize.size());
    sampleSize.push_back((double)nbs);

    canvas->cd(1);
    TH1F *scaleHist = new TH1F(Form("scale%d",nbs),";BAO Scale;Trials",50,0.75,1.25);
    scaleHist->SetLineColor(color);
    tree->Draw(Form("scale>>scale%d",nbs),"",(first ? "":"same"));
    scaleError.push_back(scaleHist->GetRMS());
    
    canvas->cd(2);
    TH1F *ampHist = new TH1F(Form("amp%d",nbs),";BAO Amplitude;Trials",50,-0.5,2.6);
    ampHist->SetLineColor(color);
    tree->Draw(Form("amp>>amp%d",nbs),"",(first ? "":"same"));
    ampError.push_back(ampHist->GetRMS());

    canvas->cd(3);
    TMarker *marker = new TMarker(nbs,scaleHist->GetRMS(),20);
    marker->SetMarkerColor(color);
    marker->Draw();
    marker = new TMarker(nbs,ampHist->GetRMS(),22);
    marker->SetMarkerColor(color);
    marker->Draw();
}

void plotBaoScaling() {
    // Initialize graphics options.
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // Initialize the drawing canvas.
    canvas = new TCanvas("canvas","canvas",1200,400);
    canvas->Divide(3,1);
    
    canvas->cd(3);
    canvas->GetPad(3)->SetLogx();
    canvas->GetPad(3)->SetLogy();
    canvas->GetPad(3)->SetGridx();
    canvas->GetPad(3)->SetGridy();
    TH2F *frame = new TH2F("frame",";Number of Plates;BAO Errors",1,75,3000,1,0.008,2.0);
    frame->Draw();

    sampleSize.clear();

    analyze(1980,kBlack);
    analyze(1000,kMagenta);
    analyze(500,kBlue);
    analyze(250,kRed);
    analyze(125,kGreen);
    
    canvas->cd(3);
    int nfit = sampleSize.size();
    TGraph *graph = new TGraph(nfit,&sampleSize[0],&scaleError[0]);
    TF1 *model = new TF1("model","[0]*pow(1000/x,0.5)",sampleSize[0],sampleSize[nfit-1]);
    model->SetLineWidth(0.1);
    graph->SetLineStyle(0);
    graph->Fit(model,"W");
    graph->Draw("P");
    model->SetParameter(0,7*model->GetParameter(0));
    model->Draw("Lsame");
}