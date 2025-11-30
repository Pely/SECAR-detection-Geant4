#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <string>

void  angles()
{
    // Open the ROOT file for reading
    // TFile* file = new TFile("../bld/cfsource_nochamb.root", "READ");
    TFile* file = new TFile("../bld_reac/Fe_pn_neutrons.root", "READ");

    // Check if the file is open
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open the ROOT file." << std::endl;
        return;
    }

    // Get the two trees from the file
    TTree* tree_det = dynamic_cast<TTree*>(file->Get("LS_reach"));
    TTree* tree_tar = dynamic_cast<TTree*>(file->Get("Target"));

    // Check if the trees exist
    if (!tree_det || !tree_tar) {
        std::cerr << "Error: One or both of the trees do not exist in the ROOT file." << std::endl;
        file->Close();
        return;
    }

    // Variables to store positions from each tree
    Double_t posX0, posY0, posZ0, posX1, posY1, posZ1, posX2, posY2, posZ2, posX3, posY3, posZ3, posX4, posY4, posZ4;
    Double_t r, theta, thetaDegrees, phi, phiDegrees;

    // Set branch addresses for tree_det
    tree_det->SetBranchAddress("posX_0", &posX1);
    tree_det->SetBranchAddress("posY_0", &posY1);
    tree_det->SetBranchAddress("posZ_0", &posZ1);

    tree_det->SetBranchAddress("posX_1", &posX2);
    tree_det->SetBranchAddress("posY_1", &posY2);
    tree_det->SetBranchAddress("posZ_1", &posZ2);

    tree_det->SetBranchAddress("posX_2", &posX3);
    tree_det->SetBranchAddress("posY_2", &posY3);
    tree_det->SetBranchAddress("posZ_2", &posZ3);

    tree_det->SetBranchAddress("posX_3", &posX4);
    tree_det->SetBranchAddress("posY_3", &posY4);
    tree_det->SetBranchAddress("posZ_3", &posZ4);

    // Set branch addresses for tree_tar
    tree_tar->SetBranchAddress("posX", &posX0);
    tree_tar->SetBranchAddress("posY", &posY0);
    tree_tar->SetBranchAddress("posZ", &posZ0);

    // Create histograms to analyze positions
    // TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
    // canvas->Divide(1,2);
    TH1F* theta_h0 = new TH1F("hist_th_0", "Theta in deg (LS0)", 100,0,50);
    TH1F* phi_h0 = new TH1F("hist_ph_0", "Phi in deg (LS0)", 720,-360,360);
    TH1F* theta_h1 = new TH1F("hist_th_1", "Theta in deg (LS1)", 100,0,50);
    TH1F* phi_h1 = new TH1F("hist_ph_1", "Phi in deg (LS1)", 720,-360,360);
    TH1F* theta_h2 = new TH1F("hist_th_2", "Theta in deg (LS2)", 100,0,50);
    TH1F* phi_h2 = new TH1F("hist_ph_2", "Phi in deg (LS2)", 720,-360,360);
    TH1F* theta_h3  = new TH1F("hist_th_3", "Theta in deg (LS3)", 100,0,50);
    TH1F* phi_h3 = new TH1F("hist_ph_3", "Phi in deg (LS3)", 720,-360,360);

    // Loop over entries in tree_det
    for (Long64_t i = 0; i < tree_det->GetEntries(); i++) 
    {
        tree_tar->GetEntry(i);
        tree_det->GetEntry(i);
        
        if(posZ0>-32&posX1>5&posY1>0&posZ1>100)
        {  
            r = TMath::Sqrt(TMath::Power((posX0-posX1),2)+TMath::Power((posY0-posY1),2)+TMath::Power((posZ0-posZ1),2));
            theta = TMath::ACos((posZ1-posZ0)/r);
            thetaDegrees = TMath::RadToDeg()*theta;
            theta_h0->Fill(thetaDegrees);
            phi = TMath::ATan2((posY1-posY0),(posX1-posX0));
            phiDegrees = TMath::RadToDeg()*phi;
            phi_h0->Fill(phiDegrees);
            //cout<<posZ0-posZ1<<" "<<TMath::ACos((posZ0-posZ1))<<endl;
        }
        if(posZ0>-32&posX2>5&posY2>0&posZ2>100)
        {  
            r = TMath::Sqrt(TMath::Power((posX0-posX2),2)+TMath::Power((posY0-posY2),2)+TMath::Power((posZ0-posZ2),2));
            theta = TMath::ACos((posZ2-posZ0)/r);
            thetaDegrees = TMath::RadToDeg()*theta;
            theta_h1->Fill(thetaDegrees);
            phi = TMath::ATan2((posY2-posY0),(posX2-posX0));
            phiDegrees = TMath::RadToDeg()*phi;
            phi_h1->Fill(phiDegrees);
        }
        if(posZ0>-32&posX3<0&posY3>0&posZ3>100)
        {  
            r = TMath::Sqrt(TMath::Power((posX0-posX3),2)+TMath::Power((posY0-posY3),2)+TMath::Power((posZ0-posZ3),2));
            theta = TMath::ACos((posZ3-posZ0)/r);
            thetaDegrees = TMath::RadToDeg()*theta;
            theta_h2->Fill(thetaDegrees);
            phi = TMath::ATan2((posY3-posY0),(posX0-posX3));
            phiDegrees = TMath::RadToDeg()*phi;
            phi_h2->Fill(phiDegrees);
        }
        if(posZ0>-32&posX4<0&posY4>0&posZ4>100)
        {  
            r = TMath::Sqrt(TMath::Power((posX0-posX4),2)+TMath::Power((posY0-posY4),2)+TMath::Power((posZ0-posZ4),2));
            theta = TMath::ACos((posZ4-posZ0)/r);
            thetaDegrees = TMath::RadToDeg()*theta;
            theta_h3->Fill(thetaDegrees);
            phi = TMath::ATan2((posY4-posY0),(posX0-posX4));
            phiDegrees = TMath::RadToDeg()*phi;
            phi_h3->Fill(phiDegrees);
        }
    }
    //canvas->Update();

    TFile *fout = new TFile("angular_coverage_pn1_LS.root", "RECREATE");
    theta_h0->Write();
    phi_h0->Write();
    theta_h1->Write();
    phi_h1->Write();
    theta_h2->Write();
    phi_h2->Write();
    theta_h3->Write();
    phi_h3->Write();
    

    // Close the ROOT file
    file->Close();
    // fout->Close();
}




