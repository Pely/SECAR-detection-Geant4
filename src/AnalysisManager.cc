
//********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Author: Pelagia Tsintari, pelagia.tsin@gmail.com
//
//

#include <stdlib.h>
#include "AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"


AnalysisManager::AnalysisManager() 
{
  factoryOn = false;
}

AnalysisManager::~AnalysisManager() 
{
  delete G4AnalysisManager::Instance();
}

void AnalysisManager::SetNeutronDetectors(G4int NeutronDetNum) { LScinNum = NeutronDetNum; }

void AnalysisManager::Book() 
{ 
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager->SetVerboseLevel(1);

  // Create a root file
  G4String fileName = "output.root";

  G4bool fileOpen = manager->OpenFile(fileName);
  if (!fileOpen) {
    G4cout << "\n---> AnalysisManager::book(): cannot open " 
           << fileName << G4endl;
    return;
  }

  manager->SetFirstNtupleId(1);

  manager -> CreateNtuple("DSSD", "DSSD variables");
  manager -> CreateNtupleDColumn(1,"totalEdep");
  manager -> CreateNtupleDColumn(1,"ekin");
  manager -> CreateNtupleDColumn(1,"A");
  manager -> CreateNtupleDColumn(1,"Z");
  manager -> CreateNtupleDColumn(1,"theta");
  manager -> CreateNtupleDColumn(1,"phi");
  manager -> CreateNtupleDColumn(1,"posx");
  manager -> CreateNtupleDColumn(1,"posy");
  manager -> FinishNtuple();

  manager -> CreateNtuple("IC", "IC variables");
  manager -> CreateNtupleDColumn(2,"totalEdep");
  manager -> CreateNtupleDColumn(2,"ekin");
  manager -> CreateNtupleDColumn(2,"A");
  manager -> CreateNtupleDColumn(2,"Z");
  manager -> FinishNtuple();

  manager -> CreateNtuple("Target", "Target variables");
  manager -> CreateNtupleDColumn(3,"Z");
  manager -> CreateNtupleDColumn(3,"A");
  manager -> CreateNtupleDColumn(3,"ekin");
  manager -> CreateNtupleDColumn(3,"edep");
  manager -> CreateNtupleDColumn(3,"posX");
  manager -> CreateNtupleDColumn(3,"posY");
  manager -> CreateNtupleDColumn(3,"posZ");
  manager -> CreateNtupleDColumn(3,"theta");
  manager -> CreateNtupleDColumn(3,"phi");
  manager -> CreateNtupleDColumn(3,"ekin_str");
  manager -> CreateNtupleDColumn(3,"posX_str");
  manager -> CreateNtupleDColumn(3,"posY_str");
  manager -> CreateNtupleDColumn(3,"posZ_str");
  manager -> CreateNtupleDColumn(3,"theta_str");
  manager -> CreateNtupleDColumn(3,"phi_str");
  manager -> CreateNtupleDColumn(3,"ToF");
  manager -> CreateNtupleDColumn(3,"ekin_n");
  manager -> CreateNtupleDColumn(3,"theta_n");
  manager -> FinishNtuple();

  manager -> CreateNtuple("LS_det", "LScinD variables");
  for(G4int i = 0; i<LScinNum; i++)
  {
    auto n = std::to_string(i);
    manager -> CreateNtupleDColumn(4,"totalEdep_"+n);
    manager -> CreateNtupleDColumn(4,"Z_"+n);
    manager -> CreateNtupleDColumn(4,"A_"+n);
  }
  manager -> FinishNtuple();

  manager -> CreateNtuple("LS_reach", "LScinR variables");
  manager -> CreateNtupleDColumn(5,"ekin_n_tar");
  manager -> CreateNtupleDColumn(5,"theta_n_tar");
  for(G4int i = 0; i<LScinNum; i++)
  {
    auto n = std::to_string(i);
    manager -> CreateNtupleDColumn(5,"ekin_"+n);
    manager -> CreateNtupleDColumn(5,"ToF_"+n);
    manager -> CreateNtupleDColumn(5,"Z_"+n);
    manager -> CreateNtupleDColumn(5,"A_"+n);
    manager -> CreateNtupleDColumn(5,"posX_"+n);
    manager -> CreateNtupleDColumn(5,"posY_"+n);
    manager -> CreateNtupleDColumn(5,"posZ_"+n);
    manager -> CreateNtupleDColumn(5,"theta_"+n);
    manager -> CreateNtupleDColumn(5,"phi_"+n);
  }
  manager -> FinishNtuple();

  manager -> CreateNtuple("Target_LS", "Target-LS coincidence");
  manager -> CreateNtupleDColumn(6,"ekin_tar");
  manager -> CreateNtupleDColumn(6,"theta_tar");
  manager -> CreateNtupleDColumn(6,"ekin_n_tar");
  manager -> CreateNtupleDColumn(6,"theta_n_tar");
  for(G4int i = 0; i<LScinNum; i++)
  {
    auto n = std::to_string(i);
    manager -> CreateNtupleDColumn(6,"ekin_"+n);
    manager -> CreateNtupleDColumn(6,"ToF_"+n);
    manager -> CreateNtupleDColumn(6,"Z_"+n);
    manager -> CreateNtupleDColumn(6,"A_"+n);
    manager -> CreateNtupleDColumn(6,"posX_"+n);
    manager -> CreateNtupleDColumn(6,"posY_"+n);
    manager -> CreateNtupleDColumn(6,"posZ_"+n);
    manager -> CreateNtupleDColumn(6,"theta_"+n);
    manager -> CreateNtupleDColumn(6,"phi_"+n);
    manager -> CreateNtupleDColumn(6,"totalEdep_"+n);
  }
  manager -> FinishNtuple();

  manager -> CreateNtuple("SiMon", "SiMon variables");
  manager -> CreateNtupleDColumn(7,"totalEdep");
  manager -> CreateNtupleDColumn(7,"ekin");
  manager -> CreateNtupleDColumn(7,"A");
  manager -> CreateNtupleDColumn(7,"Z");
  manager -> FinishNtuple();

  factoryOn = true;    
}

void AnalysisManager::Detector_DSSD(G4double edep, G4double energy,  G4double AA, G4double ZZ, G4double th, G4double phi, G4double posx, G4double posy)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager -> FillNtupleDColumn(1, 0, edep);
  manager -> FillNtupleDColumn(1, 1, energy);
  manager -> FillNtupleDColumn(1, 2, AA);
  manager -> FillNtupleDColumn(1, 3, ZZ);
  manager -> FillNtupleDColumn(1, 4, th);
  manager -> FillNtupleDColumn(1, 5, phi);
  manager -> FillNtupleDColumn(1, 6, posx);
  manager -> FillNtupleDColumn(1, 7, posy);
  manager -> AddNtupleRow(1); 
}

void AnalysisManager::Detector_IC(G4double edep, G4double energy, G4double AA, G4double ZZ)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager -> FillNtupleDColumn(2, 0, edep);
  manager -> FillNtupleDColumn(2, 1, energy);
  manager -> FillNtupleDColumn(2, 2, AA);
  manager -> FillNtupleDColumn(2, 3, ZZ);
  manager -> AddNtupleRow(2); 
}

void AnalysisManager::Target(G4int ZZ, G4int AA, G4double energy, G4double edep, G4double posx, G4double posy, G4double posz, G4double th, G4double phi, G4double energy_str, G4double posx_str, G4double posy_str, G4double posz_str, G4double th_str, G4double phi_str, G4double tof, G4double e_n, G4double th_n)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager -> FillNtupleDColumn(3, 0, ZZ);
  manager -> FillNtupleDColumn(3, 1, AA);
  manager -> FillNtupleDColumn(3, 2, energy);
  manager -> FillNtupleDColumn(3, 3, edep);
  manager -> FillNtupleDColumn(3, 4, posx);
  manager -> FillNtupleDColumn(3, 5, posy);
  manager -> FillNtupleDColumn(3, 6, posz);
  manager -> FillNtupleDColumn(3, 7, th);
  manager -> FillNtupleDColumn(3, 8, phi);
  manager -> FillNtupleDColumn(3, 9, energy_str);
  manager -> FillNtupleDColumn(3, 10, posx_str);
  manager -> FillNtupleDColumn(3, 11, posy_str);
  manager -> FillNtupleDColumn(3, 12, posz_str);
  manager -> FillNtupleDColumn(3, 13, th_str);
  manager -> FillNtupleDColumn(3, 14, phi_str);
  manager -> FillNtupleDColumn(3, 15, tof);
  manager -> FillNtupleDColumn(3, 16, e_n);
  manager -> FillNtupleDColumn(3, 17, th_n);
  manager -> AddNtupleRow(3); 
}

void AnalysisManager::Detected_LS(std::vector<G4double> edep, std::vector<G4int> z, std::vector<G4int> a)
{ 
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  for(G4int i = 0; i<LScinNum; i++)
  { 
    G4int j=3; //The number of the variables 
    manager -> FillNtupleDColumn(4,  i*j,    edep[i]);
    manager -> FillNtupleDColumn(4, (i*j)+1, z[i]);
    manager -> FillNtupleDColumn(4, (i*j)+2, a[i]);
  }
  manager -> AddNtupleRow(4);
}

void AnalysisManager::Reached_LS(std::vector<G4double> energy, std::vector<G4double> tof, std::vector<G4int> zz, std::vector<G4int> aa, std::vector<G4double> posx, std::vector<G4double> posy, std::vector<G4double> posz, std::vector<G4double> th, std::vector<G4double> phi, G4double e_nt, G4double th_nt)
{ 
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager -> FillNtupleDColumn(5, 0, e_nt);
  manager -> FillNtupleDColumn(5, 1, th_nt);
  for(G4int i = 0; i<LScinNum; i++)
  { 
    G4int j=9; //The number of the variables 
    manager -> FillNtupleDColumn(5, (i*j)+2, energy[i]);
    manager -> FillNtupleDColumn(5, (i*j)+3, tof[i]);
    manager -> FillNtupleDColumn(5, (i*j)+4, zz[i]);
    manager -> FillNtupleDColumn(5, (i*j)+5, aa[i]);
    manager -> FillNtupleDColumn(5, (i*j)+6, posx[i]);
    manager -> FillNtupleDColumn(5, (i*j)+7, posy[i]);
    manager -> FillNtupleDColumn(5, (i*j)+8, posz[i]);
    manager -> FillNtupleDColumn(5, (i*j)+9, th[i]);
    manager -> FillNtupleDColumn(5, (i*j)+10, phi[i]);
    // if(aa[i]==1&zz[i]==0)G4cout<<"Analys "<<" "<<energy[0]<<" "<<energy[1]<<" "<<energy[2]<<" "<<energy[3]<<G4endl;
  }
  manager -> AddNtupleRow(5);
}

void AnalysisManager::Target_LS_coinc(G4double energy_t, G4double th_t, G4double e_nt, G4double th_nt, std::vector<G4double> energy, std::vector<G4double> tof, std::vector<G4int> zz, std::vector<G4int> aa, std::vector<G4double> posx, std::vector<G4double> posy, std::vector<G4double> posz, std::vector<G4double> th, std::vector<G4double> phi, std::vector<G4double> edep)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager -> FillNtupleDColumn(6, 0, energy_t);
  manager -> FillNtupleDColumn(6, 1, th_t);
  manager -> FillNtupleDColumn(6, 2, e_nt);
  manager -> FillNtupleDColumn(6, 3, th_nt);
  for(G4int i = 0; i<LScinNum; i++)
  { 
    G4int j=10; //The number of the variables 
    manager -> FillNtupleDColumn(6, (i*j)+4, energy[i]);
    manager -> FillNtupleDColumn(6, (i*j)+5, tof[i]);
    manager -> FillNtupleDColumn(6, (i*j)+6, zz[i]);
    manager -> FillNtupleDColumn(6, (i*j)+7, aa[i]);
    manager -> FillNtupleDColumn(6, (i*j)+8, posx[i]);
    manager -> FillNtupleDColumn(6, (i*j)+9, posy[i]);
    manager -> FillNtupleDColumn(6, (i*j)+10, posz[i]);
    manager -> FillNtupleDColumn(6, (i*j)+11, th[i]);
    manager -> FillNtupleDColumn(6, (i*j)+12, phi[i]);
    manager -> FillNtupleDColumn(6, (i*j)+13, edep[i]);
  }
  manager -> AddNtupleRow(6);
}

void AnalysisManager::Detector_SiMon(G4double edep, G4double energy, G4double AA, G4double ZZ)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetDefaultFileType("root");
  manager -> FillNtupleDColumn(7, 0, edep);
  manager -> FillNtupleDColumn(7, 1, energy);
  manager -> FillNtupleDColumn(7, 2, AA);
  manager -> FillNtupleDColumn(7, 3, ZZ);
  manager -> AddNtupleRow(7); 
}

void AnalysisManager::Finish() 
{   
 if (factoryOn) 
   {
    G4AnalysisManager* manager = G4AnalysisManager::Instance(); 
    manager->SetDefaultFileType("root");
    manager -> Write();
    manager -> CloseFile();  

    delete G4AnalysisManager::Instance();
   }
}








