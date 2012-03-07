// Created 6-Mar-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macro to plot the scaling of BAO parameter errors with bootstrap sample size.

TCanvas *canvas, *canvas2, *canvas3;

int count = 0;
int colors[15] = { 1,2,3,4,5,6,7,8,9,46,28,30,33,38,41 };

TH1F *scaleAll, *ampAll;

void analyze(int index, const char *pattern) {
    
    int color = colors[count];

    // Load the tree
    TFile *file = new TFile(Form(pattern,index),"READ");
    TTree *tree = (TTree*)file->Get("bs");

    // Plot a scatter plot of each bootstrap sample in (scale,amp)
    canvas->cd(1);
    TH2F *scaleAmpHist = new TH2F(Form("scaleamp%d",count),";BAO Scale;BAO Amplitude",
        100,0.65,1.35,100,-1.,5.);
    scaleAmpHist->SetMarkerColor(color);
    scaleAmpHist->SetMarkerStyle(4);
    scaleAmpHist->SetMarkerSize(0.4);
    tree->Draw(Form("amp:scale>>scaleamp%d",count),"","goff");
    scaleAmpHist->Draw((count ? "same":""));
    
    // Plot a histogram of bootstrap scale values (not plotted)
    TH1F *scaleHist = new TH1F(Form("scale%d",count),";BAO Scale;Bootstrap Trials",50,0.65,1.35);
    tree->Draw(Form("scale>>scale%d",count),"","goff");
    
    // Plot a histogram of bootstrap amplitude values (not plotted)
    TH1F *ampHist = new TH1F(Form("amp%d",count),";BAO Amplitude;Bootstrap Trials",50,-1,5);
    tree->Draw(Form("amp>>amp%d",count),"","goff");
    
    // Plot per-realization bootstrap distributions of scale and amplitude.
    canvas2->cd(count+1);
    scaleHist->SetFillColor(color);
    scaleHist->Draw();
    canvas3->cd(count+1);
    ampHist->SetFillColor(color);
    ampHist->Draw();
    
    // Accumulate histograms of bootstrap scale and amplitude for plotting at the end.
    if(0 == count) {
        scaleAll = (TH1F*)scaleHist->Clone("scaleAll");
        ampAll = (TH1F*)ampHist->Clone("ampAll");
    }
    else {
        scaleAll->Add(scaleHist);
        ampAll->Add(ampHist);
    }
    
    // Plot a marker with error bars in (scale,amp) showing the mean and sigma of this sample.
    canvas->cd(4);
    double ampMean[1],scaleMean[1],ampRms[1],scaleRms[1];
    ampMean[0] = ampHist->GetMean();
    scaleMean[0] = scaleHist->GetMean();
    ampRms[0] = ampHist->GetRMS();
    scaleRms[0] = scaleHist->GetRMS();
    TGraphErrors *graph = new TGraphErrors(1,scaleMean,ampMean,scaleRms,ampRms);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(4);
    graph->SetLineColor(color);
    graph->Draw("P");
    
    count++;
}

void plotBaoSamples() {
    // Initialize graphics options.
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // Initialize the drawing canvas.
    canvas = new TCanvas("canvas","canvas",1000,1000);
    canvas->Divide(2,2);
    
    canvas2 = new TCanvas("canvas2","canvas2",1500,900);
    canvas2->Divide(5,3,0,0);
    canvas3 = new TCanvas("canvas3","canvas3",1500,900);
    canvas3->Divide(5,3,0,0);
    
    canvas->cd(4);
    TH2F *frame4 = new TH2F("frame4",";BAO Scale;BAO Amplitude",1,0.65,1.35,1,-1.,5.);
    frame4->Draw();
    
    count = 0;    
    const char *no_noise = "results/delta_diag_%d.root";
    const char *with_noise = "results/delta_diag_n_%d.root";

    for(int index = 1; index <= 15; ++index) {
        analyze(index, no_noise);
    }
    
    canvas->cd(2);
    ampAll->SetLineWidth(1.5);
    ampAll->SetFillColor(kBlue-8);
    ampAll->Draw();
    TText *ampLabel = new TLatex(0.55,0.8,Form("%.3f #pm %.3f",ampAll->GetMean(),ampAll->GetRMS()));
    ampLabel->SetTextColor(kBlue-8);
    ampLabel->SetNDC();
    ampLabel->Draw();

    canvas->cd(3);
    scaleAll->SetLineWidth(1.5);
    scaleAll->SetFillColor(kBlue-8);
    scaleAll->Draw();    
    TText *scaleLabel = new TLatex(0.15,0.8,Form("%.3f #pm %.3f",scaleAll->GetMean(),scaleAll->GetRMS()));
    scaleLabel->SetTextColor(kBlue-8);
    scaleLabel->SetNDC();
    scaleLabel->Draw();
}
