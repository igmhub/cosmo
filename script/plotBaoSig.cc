// Created 13-Mar-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// ROOT macro to plot the bootstrap significance estimates for 15 data realizations.

TCanvas *canvas;

int count = 0;
int colors[15] = { 1,2,3,4,5,6,7,8,9,46,28,30,33,38,41 };

void analyze(int index, const char *pattern1, const char *pattern2) {
    
    int color1 = colors[count], color2 = (0==count)? kRed:kBlack;
    
    // Load the trees.
    TFile *file1 = new TFile(Form(pattern1,index),"READ");
    TTree *tree1 = (TTree*)file1->Get("bs");

    TFile *file2 = new TFile(Form(pattern2,index),"READ");
    TTree *tree2 = (TTree*)file2->Get("bs");

    // Plot per-realization bootstrap distributions of p1 and p2.
    int pad = count+1;
    canvas->cd(pad);

    // Create empty histograms.
    int nbins = 50;
    TH1F *hist1 = new TH1F(Form("hist1_%d",count),";BAO Amplitude;",nbins,-5.,5.);
    TH1F *hist2 = (TH1F*)hist1->Clone(Form("hist2_%d",count));

    // Draw the first histogram with a filled background.
    tree1->Draw(Form("amp>>hist1_%d",count));
    hist1->SetFillColor(color1);

    // Draw the second histogram with a hatched background.
    tree2->Draw(Form("amp>>hist2_%d",count),"","same");
    hist2->SetLineColor(color2);
    hist2->SetFillColor(color2);
    hist2->SetFillStyle(3004);
    
    // Adjust the vertical range for the second histogram, if necessary.
    double ymax = hist2->GetMaximum();
    if(ymax > hist1->GetMaximum()) {
        hist1->SetMaximum(1.05*ymax);
    }
    else {
        ymax = hist1->GetMaximum();
    }

    // Fetch the integrals of the histograms.
    double *integral1 = hist1->GetIntegral();
    double *integral2 = hist2->GetIntegral();

    // Fetch bin centers and calculate the integrated confidence level.
    double sumCL = 0;
    double *binCenter = new double[nbins];
    for(int bin = 0; bin < nbins; ++bin) {
        binCenter[bin] = hist2->GetXaxis()->GetBinCenter(bin+1);
        double binCL = integral2[bin];
        integral1[bin] = ymax*(1-integral1[bin]);
        integral2[bin] = ymax*binCL;
        sumCL += binCL*hist1->GetBinContent(bin+1);
    }

    // Draw the integral as a graph.
    TGraph *graph2 = new TGraph(nbins,binCenter,integral2);
    graph2->SetLineColor(color2);
    graph2->SetLineWidth(1.5);
    graph2->Draw("L");
    TGraph *graph1 = new TGraph(nbins,binCenter,integral1);
    graph1->SetLineColor(color2);
    graph1->SetLineWidth(1.5);
    graph1->SetLineStyle(kDashed);
    graph1->Draw("L");
    
    // Interpolate the hist2 integral at the hist1 mean.
    double hist1Mean = hist1->GetMean();
    double meanCL = graph2->Eval(hist1Mean)/ymax;
    double zeroFrac = graph1->Eval(0)/ymax;
    
    // Calculate the average CL
    double avgCL = sumCL/hist1->GetEntries();
    std::cout << index << ' ' << avgCL << ' ' << meanCL << std::endl;
    
    // Draw the hist1 mean.
    canvas->Update();
    ymax = canvas->GetPad(pad)->GetUymax();
    TLine *line = new TLine(hist1Mean,0.,hist1Mean,ymax);
    line->SetLineWidth(1.5);
    line->SetLineColor(color2);
    line->Draw();
    TLine *line2 = new TLine(0,0.,0,ymax);
    line2->SetLineWidth(1.5);
    line2->SetLineColor(color2);
    line2->SetLineStyle(kDashed);
    line2->Draw();

    // Draw the CL in the top-left corner.
    double xLabel = (1 == (pad % 5)) ? 0.15 : 0.06;
    TText *label2 = new TLatex(xLabel,0.7,Form("%.1f%%",100*zeroFrac));
    label2->SetTextFont(33);
    label2->SetTextAlign(13);
    label2->SetTextSize(24);
    label2->SetNDC();
    label2->Draw();
    TText *label1 = new TLatex(0.95,0.7,Form("%.1f%%",100*meanCL));
    label1->SetTextFont(33);
    label1->SetTextAlign(33);
    label1->SetTextSize(24);
    label1->SetNDC();
    label1->Draw();
    
    count++;
}

/**
Patterns are a printf string with a single %d field that will be substituted with 1-15 to
read results for each realization, e.g.

  results/delta_diag_%d.root
  results/delta_diag_n_%d.root
  results/fix_bao_no_noise_%d.root

**/
void plotBaoSig(const char *pattern, const char *nullPattern) {
    // Initialize graphics options.
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // Initialize the drawing canvases.
    canvas = new TCanvas("canvas","canvas",1500,900);
    canvas->Divide(5,3,0,0);
    
    count = 0;
    for(int index = 1; index <= 15; ++index) {
        analyze(index, pattern, nullPattern);
    }
}
