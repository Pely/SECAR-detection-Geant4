#include <TH2F.h>
#include <fstream>

void HistEntriesToTxt() {

    // Open the ROOT file
    TFile* file = new TFile("../bld_reac/cobalt_ch2_withchamb.root", "READ");

    // Check if the file is open
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open the ROOT file." << std::endl;
        return;
    }

    // Get the Target tree from the file
    TTree* tree = dynamic_cast<TTree*>(file->Get("Target"));
    TTree* tree_dssd = dynamic_cast<TTree*>(file->Get("DSSD"));
    tree->AddFriend(tree_dssd);

    TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
    TH2F* hist2D = new TH2F("hist2D", "2D Histogram", 300,0,30,20,192,202);
    //TH1D* hist1D = new TH1D("hist1D", "1D Histogram", 100, 0, 10);
    // tree->Draw("Ekin:theta*TMath::RadToDeg()>>hist2D", "(Z==0 && A==1) & (DSSD.Z==27 & DSSD.A==58) & theta<0.5", "colz");
    tree->Draw("Ekin:theta*1000>>hist2D", "", "colz");
    //tree->Draw("theta>>hist1D", "Z==0 && A==1 & theta<=0.44 & Ekin<8 & Ekin>1");
    // tree->Draw("Ekin>>hist1D", "Z==0 && A==1 & theta<=0.44 & Ekin<10 & Ekin>1");

    // Create a text file for saving entries
    std::ofstream outputFile("Co_entries.txt");

    // Loop over the bins and save the bin contents (energy, angle, and count)
    for (int ybin = 1; ybin <= hist2D->GetNbinsY(); ++ybin) {
        for (int xbin = 1; xbin <= hist2D->GetNbinsX(); ++xbin) {
            double angle = hist2D->GetXaxis()->GetBinCenter(xbin);
            double energy = hist2D->GetYaxis()->GetBinCenter(ybin);
            int count = static_cast<int>(hist2D->GetBinContent(xbin, ybin));
            // int count = static_cast<int>(hist1D->GetBinContent(xbin));

            // Write data to the text file
            if(count>0) outputFile << energy << "\t" << angle << std::endl; //<< "\t" << count 
        }
    }

    canvas->Draw();
    canvas->Update();

    // Close the files
    outputFile.close();
    //file->Close();
}