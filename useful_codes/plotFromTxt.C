#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>

void plotFromTxt() {
    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);

    // Create a histogram
    TH2D* hist = new TH2D("hist", "Data from Txt File", 250, 0, 25, 100, 0, 10);

    // Open the text file
    std::ifstream inputFile("../neutronInput/Fe_pn_entries.txt");

    // Check if the file is open
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open the data file." << std::endl;
        return;
    }

    // Read data from the file and fill the histogram
    double x, y, z;
    while (inputFile >> x >> y >> z) {
        hist->Fill(y,x);
    }

    // Close the text file
    inputFile.close();

    // Draw the histogram
    hist->Draw("colz");

    // Update the canvas
    canvas->Update();
}
