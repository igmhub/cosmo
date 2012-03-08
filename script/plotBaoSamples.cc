// Created 6-Mar-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macro to plot the scaling of BAO parameter errors with bootstrap sample size.

TCanvas *canvas, *canvas2, *canvas3;

int count = 0;
int colors[15] = { 1,2,3,4,5,6,7,8,9,46,28,30,33,38,41 };

char name[9][16] = { "alpha", "bias", "beta", "amp", "scale", "xio", "a0", "a1", "a2" };

int bins[9]   = {   50,    50,  50,  50,   50,    50,  50,  50,  50   };
double xlo[9] = { -0.5, -0.05, -1., -1., 0.65, -1e-3, -3., -10., -10. };
double xhi[9] = {  9.5,  0.50, 10.,  5., 1.35,  1e-3,  9.,  10.,  10. };

TH1F *scaleAll, *ampAll;

void analyze(int index, const char *pattern, int p1, int p2) {
    
    int color = colors[count];
    
    // Load the tree
    TFile *file = new TFile(Form(pattern,index),"READ");
    TTree *tree = (TTree*)file->Get("bs");

    // Plot a scatter plot of each bootstrap sample in (scale,amp)
    canvas->cd(1);
    TH2F *scaleAmpHist = new TH2F(Form("%s%s%d",name[p1],name[p2],count),Form(";%s;%s",name[p1],name[p2]),
        100,xlo[p1],xhi[p1],100,xlo[p2],xhi[p2]);
    scaleAmpHist->SetMarkerColor(color);
    scaleAmpHist->SetMarkerStyle(4);
    scaleAmpHist->SetMarkerSize(0.4);
    tree->Draw(Form("%s:%s>>%s%s%d",name[p2],name[p1],name[p1],name[p2],count),"","goff");
    scaleAmpHist->Draw((count ? "same":""));
    
    // Plot a histogram of bootstrap scale values (not plotted)
    TH1F *scaleHist = new TH1F(Form("%s%d",name[p1],count),Form(";%s;Bootstrap Trials",name[p1]),
        bins[p1],xlo[p1],xhi[p1]);
    tree->Draw(Form("%s>>%s%d",name[p1],name[p1],count),"","goff");
    
    // Plot a histogram of bootstrap amplitude values (not plotted)
    TH1F *ampHist = new TH1F(Form("%s%d",name[p2],count),Form(";%s;Bootstrap Trials",name[p2]),
        bins[p2],xlo[p2],xhi[p2]);
    tree->Draw(Form("%s>>%s%d",name[p2],name[p2],count),"","goff");
    
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

/**
Pattern is a printf string with a single %d field that will be substituted with 1-15 to
read results for each realization, e.g.

  results/delta_diag_%d.root
  results/delta_diag_n_%d.root
  results/fix_bao_no_noise_%d.root

Select parameter indices p1,p2 according to:

0   ||     Alpha
1   ||      Bias
2   ||      Beta
3   ||  BAO Ampl
4   || BAO Scale
5   ||    BB xio
6   ||     BB a0
7   ||     BB a1
8   ||     BB a2

**/
void plotBaoSamples(const char *pattern, int p1 = 4, int p2 = 3) {
    // Initialize graphics options.
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // Initialize the drawing canvases.
    canvas = new TCanvas("canvas","canvas",1000,1000);
    canvas->Divide(2,2);
    canvas2 = new TCanvas("canvas2","canvas2",1500,900);
    canvas2->Divide(5,3,0,0);
    canvas3 = new TCanvas("canvas3","canvas3",1500,900);
    canvas3->Divide(5,3,0,0);
    
    canvas->cd(4);
    TH2F *frame4 = new TH2F("frame4",Form(";%s;%s",name[p1],name[p2]),
        1,xlo[p1],xhi[p1],1,xlo[p2],xhi[p2]);
    frame4->Draw();
    
    count = 0;    
    const char *no_noise = "results/delta_diag_%d.root";
    const char *with_noise = "results/delta_diag_n_%d.root";

    for(int index = 1; index <= 15; ++index) {
        analyze(index, "results/fix_bao_no_noise_%d.root", p1, p2);
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
