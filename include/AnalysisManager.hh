//
// ********************************************************************
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
 

#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH 

#include "globals.hh"
#include "G4AnalysisManager.hh"

class AnalysisManager
{ 

public:
   AnalysisManager();
  ~AnalysisManager();
  
  void Book(); // booking the ROOT file

  void Target(G4int ZZ, G4int AA, G4double energy, G4double edep, G4double posx, G4double posy, G4double posz, G4double th, G4double phi, G4double energy_str, G4double posx_str, G4double posy_str, G4double posz_str, G4double th_str, G4double phi_str, G4double tof, G4double e_n, G4double th_n);

  void Detector_SiMon(G4double edep, G4double energy, G4double AA, G4double ZZ);
  
  void Detector_DSSD(G4double edep, G4double energy, G4double AA, G4double ZZ, G4double th, G4double phi, G4double posx, G4double posy);//, G4int coinc, G4double energy_k, G4double energy_n, G4double pn_x, G4double pn_y, G4double th_n);
  
  void Detector_IC(G4double edep, G4double energy, G4double AA, G4double ZZ);

  void Detected_LS(std::vector<G4double> Edep, std::vector<G4int> z, std::vector<G4int> a); 
  
  void Reached_LS(std::vector<G4double> Energy, std::vector<G4double> tof, std::vector<G4int> zz, std::vector<G4int> aa, std::vector<G4double> posx, std::vector<G4double> posy, std::vector<G4double> posz, std::vector<G4double> th, std::vector<G4double> phi, G4double e_nt, G4double th_nt); 
  
  void SetNeutronDetectors(G4int NeutronDetNum);

  void Target_LS_coinc(G4double Energy_t, G4double th_t, G4double e_nt, G4double th_nt, std::vector<G4double> Energy, std::vector<G4double> tof, std::vector<G4int> zz, std::vector<G4int> aa, std::vector<G4double> posx, std::vector<G4double> posy, std::vector<G4double> posz, std::vector<G4double> th, std::vector<G4double> phi, std::vector<G4double> Edep);
  
  void Finish();
  // Close the ROOT file with all the results stored in nutples 

private:

  G4bool factoryOn;
  G4int LScinNum;

};

#endif




