#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

void plotAndConvertToTxt() {
    auto inputFile = "../bld_reac/tarttest.root";
    //auto inputFile = "../bld_reac/neut_rec_htar_LS.root";
    auto outputFile = "../bld_reac/recoil_neut_coinc.dat";
    auto outputFile_eff = "../bld_reac/efficiency.dat";

    // Open the ROOT file
    TFile* file = new TFile(inputFile);

    // Check if the file is open successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open the input file." << std::endl;
        return;
    }

    // Get the TTree from the file
    TTree* tree = dynamic_cast<TTree*>(file->Get("Target_LS"));
    TTree* tree_t = dynamic_cast<TTree*>(file->Get("Target"));
    TTree* tree_det = dynamic_cast<TTree*>(file->Get("LS_det"));

    // Check if the tree is found
    if (!tree) {
        std::cerr << "Error: Unable to find the TTree in the ROOT file." << std::endl;
        file->Close();
        return;
    }

    TFile* outputHistFile = new TFile("../bld_reac/recoil_neut_coinc_hists_test.root", "RECREATE");

    // Create 2D histograms for the four branches
    TH2F* hist1_r = new TH2F("hist1_r", "Recoil distribution in coincidence with detected neutrons from LS0", 200,0,0.02,1500,190,205);
    TH2F* hist2_r = new TH2F("hist2_r", "Recoil distribution in coincidence with detected neutrons from LS1", 200,0,0.02,1500,190,205);
    TH2F* hist3_r = new TH2F("hist3_r", "Recoil distribution in coincidence with detected neutrons from LS2", 200,0,0.02,1500,190,205);
    TH2F* hist4_r = new TH2F("hist4_r", "Recoil distribution in coincidence with detected neutrons from LS3", 200,0,0.02,1500,190,205);

    TH2F* hist_r = new TH2F("hist_r", "Recoil distribution", 200,0,0.02,1500,190,205);
    TH2F* hist_r_n = new TH2F("hist_r_n", "Recoil distribution in coincidence with all detected neutrons", 200,0,0.02,1500,190,205);

    TH2F* hist1_n = new TH2F("hist1_n", "Neutron distribution in coincidence with detected neutrons from LS0", 500,0,0.5,800,0,8);
    TH2F* hist2_n = new TH2F("hist2_n", "Neutron distribution in coincidence with detected neutrons from LS1", 500,0,0.5,800,0,8);
    TH2F* hist3_n = new TH2F("hist3_n", "Neutron distribution in coincidence with detected neutrons from LS2", 500,0,0.5,800,0,8);
    TH2F* hist4_n = new TH2F("hist4_n", "Neutron distribution in coincidence with detected neutrons from LS3", 500,0,0.5,800,0,8);

    TH2F* hist_n = new TH2F("hist_n", "Neutron distribution", 500,0,0.5,800,0,8);
    TH2F* hist_n_n = new TH2F("hist_n_n", "Neutron distribution in coincidence with all detected neutrons", 500,0,0.5,800,0,8);

    TH1F* hist1_n_e = new TH1F("hist1_n_e", "Neutron energy in coincidence with detected neutrons from LS0", 50, 0, 10);
    TH1F* hist2_n_e = new TH1F("hist2_n_e", "Neutron energy in coincidence with detected neutrons from LS1", 50, 0, 10);
    TH1F* hist3_n_e = new TH1F("hist3_n_e", "Neutron energy in coincidence with detected neutrons from LS2", 50, 0, 10);
    TH1F* hist4_n_e = new TH1F("hist4_n_e", "Neutron energy in coincidence with detected neutrons from LS3", 50, 0, 10);

    TH1F* hist1_n_edep = new TH1F("hist1_n_edep", "Energy deposited for detected neutrons from LS0", 1000, 0, 10);
    TH1F* hist2_n_edep = new TH1F("hist2_n_edep", "Energy deposited for detected neutrons from LS1", 1000, 0, 10);
    TH1F* hist3_n_edep = new TH1F("hist3_n_edep", "Energy deposited for detected neutrons from LS2", 1000, 0, 10);
    TH1F* hist4_n_edep = new TH1F("hist4_n_edep", "Energy deposited for detected neutrons from LS3", 1000, 0, 10);

    TH1F* hist_n_e_n = new TH1F("hist_n_e_n", "Neutron energy in coincidence with all detected neutrons", 50, 0, 10);
    TH1F* hist_n_e = new TH1F("hist_n_e", "Neutron energy", 50, 0, 10);

    TGraph* scatterPlot = new TGraph();

    // Recoil histograms in coincidence with neutron detection
    tree->Draw("ekin_tar:theta_tar>>hist1_r", "((totalEdep_0>0.240&Z_0==1&A_0==1)||(totalEdep_0>3&Z_0>1&A_0>1))&&(theta_tar>0&ekin_tar>0&totalEdep_1==0&totalEdep_2==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_tar:theta_tar>>hist2_r", "((totalEdep_1>0.300&Z_1==1&A_1==1)||(totalEdep_1>3&Z_1>1&A_1>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_2==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_tar:theta_tar>>hist3_r", "((totalEdep_2>0.240&Z_2==1&A_2==1)||(totalEdep_2>3&Z_2>1&A_2>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_1==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_tar:theta_tar>>hist4_r", "((totalEdep_3>0.314&Z_3==1&A_3==1)||(totalEdep_3>3&Z_3>1&A_3>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_1==0&totalEdep_2==0)", "colz");

    tree->Draw("ekin_n_tar:theta_n_tar>>hist1_n", "((totalEdep_0>0.240&Z_0==1&A_0==1)||(totalEdep_0>3&Z_0>1&A_0>1))&&(theta_tar>0&ekin_tar>0&totalEdep_1==0&totalEdep_2==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_n_tar:theta_n_tar>>hist2_n", "((totalEdep_1>0.300&Z_1==1&A_1==1)||(totalEdep_1>3&Z_1>1&A_1>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_2==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_n_tar:theta_n_tar>>hist3_n", "((totalEdep_2>0.240&Z_2==1&A_2==1)||(totalEdep_2>3&Z_2>1&A_2>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_1==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_n_tar:theta_n_tar>>hist4_n", "((totalEdep_3>0.314&Z_3==1&A_3==1)||(totalEdep_3>3&Z_3>1&A_3>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_1==0&totalEdep_2==0)", "colz");

    // All recoils emmitted from target
    tree_t->Draw("ekin:theta>>hist_r", "", "colz"); 

    // All neutrons emmitted from target 
    tree_t->Draw("ekin_n:theta_n>>hist_n", "", "colz"); 

    // The cumulative recoils emmitted from target in coincidence with detected neutrons
    for (int i = 1; i <= hist_r_n->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist_r_n->GetNbinsY(); ++j) {
            hist_r_n->SetBinContent(i, j,
                hist1_r->GetBinContent(i, j) +
                hist2_r->GetBinContent(i, j) +
                hist3_r->GetBinContent(i, j) +
                hist4_r->GetBinContent(i, j));
        }
    }
    hist_r_n->Draw("colz");

    // The cumulative neutrons emmitted from target gated to the detected ones
    for (int i = 1; i <= hist_n_n->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist_n_n->GetNbinsY(); ++j) {
            hist_n_n->SetBinContent(i, j,
                hist1_n->GetBinContent(i, j) +
                hist2_n->GetBinContent(i, j) +
                hist3_n->GetBinContent(i, j) +
                hist4_n->GetBinContent(i, j));
        }
    }
    hist_n_n->Draw("colz");

    tree->Draw("ekin_n_tar>>hist1_n_e", "((totalEdep_0>0.240&Z_0==1&A_0==1)||(totalEdep_0>3&Z_0>1&A_0>1))&&(theta_tar>0&ekin_tar>0&totalEdep_1==0&totalEdep_2==0&totalEdep_3==0)&(ekin_1==0&ekin_2==0&ekin_3==0)", "colz");
    tree->Draw("ekin_n_tar>>hist2_n_e", "((totalEdep_1>0.300&Z_1==1&A_1==1)||(totalEdep_1>3&Z_1>1&A_1>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_2==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_n_tar>>hist3_n_e", "((totalEdep_2>0.240&Z_2==1&A_2==1)||(totalEdep_2>3&Z_2>1&A_2>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_1==0&totalEdep_3==0)", "colz");
    tree->Draw("ekin_n_tar>>hist4_n_e", "((totalEdep_3>0.314&Z_3==1&A_3==1)||(totalEdep_3>3&Z_3>1&A_3>1))&&(theta_tar>0&ekin_tar>0&totalEdep_0==0&totalEdep_1==0&totalEdep_2==0)", "colz");

    tree_det->Draw("totalEdep_0>>hist1_n_edep", "((totalEdep_0>0.240&Z_0==1&A_0==1)||(totalEdep_0>3&Z_0>1&A_0>1))&&(totalEdep_1==0&totalEdep_2==0&totalEdep_3==0)", "colz");
    tree_det->Draw("totalEdep_1>>hist2_n_edep", "((totalEdep_1>0.300&Z_1==1&A_1==1)||(totalEdep_1>3&Z_1>1&A_1>1))&&(totalEdep_0==0&totalEdep_2==0&totalEdep_3==0)", "colz");
    tree_det->Draw("totalEdep_2>>hist3_n_edep", "((totalEdep_2>0.240&Z_2==1&A_2==1)||(totalEdep_2>3&Z_2>1&A_2>1))&&(totalEdep_0==0&totalEdep_1==0&totalEdep_3==0)", "colz");
    tree_det->Draw("totalEdep_3>>hist4_n_edep", "((totalEdep_3>0.314&Z_3==1&A_3==1)||(totalEdep_3>3&Z_3>1&A_3>1))&&(totalEdep_0==0&totalEdep_1==0&totalEdep_2==0)", "colz");

    // All neutrons emmitted from target 
    tree_t->Draw("ekin_n>>hist_n_e", "", "colz");

    // The cumulative neutrons emmitted from target gated to the detected ones
    for (int i = 1; i <= hist_n_e_n->GetNbinsX(); ++i) {
        hist_n_e_n->SetBinContent(i,
        hist1_n_e->GetBinContent(i) +
        hist2_n_e->GetBinContent(i) +
        hist3_n_e->GetBinContent(i) +
        hist4_n_e->GetBinContent(i));
    }
    hist_n_e_n->Draw("colz");

    // Save histogram values to a single data file
    //std::ofstream txtFile_eff(outputFile_eff);
    for (int i = 0; i <= hist_n_e->GetNbinsX(); i++) {

        int count_d = static_cast<int>(hist_n_e_n->GetBinContent(i));
        int count_t = static_cast<int>(hist_n_e->GetBinContent(i));
        double ratio = static_cast<double>(count_d) / count_t;
        cout<<count_d<<"\t"<<count_t<<"\t"<<ratio<<endl;
        if (ratio>0 || !std::isnan(ratio)) {
            double xValue = hist_n_e->GetXaxis()->GetBinCenter(i);
            double yValue = ratio*100;
            scatterPlot->SetPoint(scatterPlot->GetN(), xValue, yValue);
            //txtFile_eff << xValue << "\t" << count_d << "\t" << count_t << "\t" << yValue << std::endl;
        }
    }
    scatterPlot->SetTitle("Total Neutron Detection Efficiency;Neutron Energy (MeV);Efficiency (%)");
    scatterPlot->SetMarkerStyle(20);
    scatterPlot->SetMarkerSize(0.5);
    scatterPlot->Draw("AP"); 
    scatterPlot->Write("LS_Efficiency");

    // Save histograms to a ROOT file
    hist1_r->Write();
    hist2_r->Write();
    hist3_r->Write();
    hist4_r->Write();
    hist_r->Write();
    hist_r_n->Write();

    hist1_n->Write();
    hist1_n_e->Write();
    hist2_n->Write();
    hist2_n_e->Write();
    hist3_n->Write();
    hist3_n_e->Write();
    hist4_n->Write();
    hist4_n_e->Write();

    hist_n->Write();
    hist_n_e->Write();
    hist_n_n->Write();
    hist_n_e_n->Write();

    hist1_n_edep->Write();
    hist2_n_edep->Write();
    hist3_n_edep->Write();
    hist4_n_edep->Write();

    // Save histogram values to a single data file
    /*std::ofstream txtFile(outputFile);

    for (int i = 1; i <= hist1_r->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist1_r->GetNbinsY(); ++j) {
            double x = hist1_r->GetXaxis()->GetBinCenter(i);
            double y = hist1_r->GetYaxis()->GetBinCenter(j);
            int count1 = static_cast<int>(hist1_r->GetBinContent(i, j));
            
            if(count1>0) txtFile << y << "\t" << x*1000 << "\t" << count1 << std::endl;
        }
    }

    for (int i = 1; i <= hist2_r->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist2_r->GetNbinsY(); ++j) {
            double x = hist2_r->GetXaxis()->GetBinCenter(i);
            double y = hist2_r->GetYaxis()->GetBinCenter(j);
            int count2 = static_cast<int>(hist2_r->GetBinContent(i, j));
            
            if(count2>0) txtFile << y << "\t" << x*1000 << "\t" << count2 << std::endl;
        }
    }
    for (int i = 1; i <= hist3_r->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist3_r->GetNbinsY(); ++j) {
            double x = hist3_r->GetXaxis()->GetBinCenter(i);
            double y = hist3_r->GetYaxis()->GetBinCenter(j);
            int count3 = static_cast<int>(hist3_r->GetBinContent(i, j));
            
            if(count3>0) txtFile << y << "\t" << x*1000 << "\t" << count3 << std::endl;
        }
    }
    for (int i = 1; i <= hist4_r->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist4_r->GetNbinsY(); ++j) {
            double x = hist4_r->GetXaxis()->GetBinCenter(i);
            double y = hist4_r->GetYaxis()->GetBinCenter(j);
            int count4 = static_cast<int>(hist4_r->GetBinContent(i, j));
            
            if(count4>0) txtFile << y << "\t" << x*1000 << "\t" << count4 << std::endl;
        }
    }
    txtFile.close();
    txtFile_eff.close();*/


    outputHistFile->Close();
    
    if (file) {
    file->Close();
    delete file;
    }
}


