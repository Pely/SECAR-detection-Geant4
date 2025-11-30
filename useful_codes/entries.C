#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <string>


void entries()
{
    char file_name[200], gates_tar[200];

    for (int j=100; j<=1400; j+=100)
    {   
        
        sprintf(file_name,"../bld/neutrons/4pi/n_%dkeV.root",j);
        TFile* file = new TFile(file_name, "READ");

        // Get the main tree (LS_det) from the file
        TTree* mainTree = dynamic_cast<TTree*>(file->Get("LS_det"));
        TTree* tarTree = dynamic_cast<TTree*>(file->Get("Target"));

        TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
   
            TH1F* hist = new TH1F("hist", "Total Edep", 160, 0, 16);
            mainTree->Draw("totalEdep_3>>hist", "(A_3>=12&&Z_3>=6&&totalEdep_3>3)||(A_3==1&&Z_3==1&&totalEdep_3>0.275)");
            int numEntries = hist->GetEntries();
            // std::cout << numEntries << std::endl;
            hist->Delete();
        

            // sprintf(gates_tar,"A==1 && Z==0 && Ekin>0.05");

            TH1F* histtar = new TH1F("histtar", "Neutrons from Target", 160, 0, 16);
            reachTree->Draw("Ekin_3>>histtar", "A_3==1 && Z_3==0");
            int numEntries_t = histtar->GetEntries();
            // std::cout << j*0.001 <<"\t"<< numEntries_t << std::endl;
            histtar->Delete();
        

        canvas->Close();
        file->Close();
    }
    for (int j=1500; j<=10000; j+=500)
    {
        // Open the ROOT file that contains your trees
        sprintf(file_name,"../bld/neutrons/4pi/n_%dkeV.root",j);
        TFile* file = new TFile(file_name, "READ");

        // Get the main tree (LS_det) from the file
        TTree* mainTree = dynamic_cast<TTree*>(file->Get("LS_det"));
        TTree* tarTree = dynamic_cast<TTree*>(file->Get("Target"));

        TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);

        TH1F* hist = new TH1F("hist", "Total Edep", 160, 0, 16);
        mainTree->Draw("totalEdep_3>>hist", "(A_3>=12&&Z_3>=6&&totalEdep_3>3)||(A_3==1&&Z_3==1&&totalEdep_3>0.275)");
        int numEntries = hist->GetEntries();
        // std::cout << numEntries << std::endl;
        hist->Delete();
       
        TH1F* histtar = new TH1F("histtar", "Neutrons from Target", 160, 0, 16);
        reachTree->Draw("Ekin_3>>histtar", "A_3==1 && Z_3==0");
        int numEntries_t = histtar->GetEntries();
        // std::cout << j << "   det   " << numEntries << "   tar   " << numEntries_t << std::endl;
        // std::cout << j*0.001 <<"\t"<< numEntries_t << std::endl;
        histtar->Delete();
 
        //}
        canvas->Close();
        file->Close();
    }
}
