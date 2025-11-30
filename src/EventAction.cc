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
//

#include "EventAction.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include <vector>
#include <fstream> 
#include <iostream>

EventAction::EventAction(AnalysisManager* analysisMan, DetectorConstruction* det)
  :G4UserEventAction(),
   detector(det)
{
  analysis = analysisMan;
  LScinNum = detector->GetLScinNum();
  LScin = detector->LScinOn();

  if(LScin)
  {
    for(G4int i = 0; i<LScinNum; i++)
    {
      LSEkin.push_back(0.0); 
      LSToF.push_back(0.0);
      ZZ.push_back(0);
      AA.push_back(0); 
      LSPosX.push_back(0.0);
      LSPosY.push_back(0.0); 
      LSPosZ.push_back(0.0); 
      LSTheta.push_back(0.0); 
      LSPhi.push_back(0.0); 

      LSEdep.push_back(0.0);
      Z.push_back(0);
      A.push_back(0);
    }
  }
  
  TarZ = 0;
  TarA = 0;
  TarEnergy = 0.0;
  TarEdep = 0.0;
  TarPosX = 0.0;
  TarPosY = 0.0; 
  TarPosZ = 0.0;
  TarTheta = 0.0; 
  TarPhi = 0.0;
  TarToF = 0.0;
  Tar_n_e = 0.0;
  Tar_n_theta = 0.0;

  StripEnergy = 0.0;
  StripPosX = 0.0; 
  StripPosY = 0.0; 
  StripPosZ = 0.0; 
  StripTheta = 0.0; 
  StripPhi = 0.0; 
}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event* evt)
{ 
  if(LScin)
  {  
    for(G4int i = 0; i<LScinNum; i++)
    {
      std::fill(LSEkin.begin(),LSEkin.end(),0.0);  
      std::fill(LSToF.begin(),LSToF.end(),0.0); 
      std::fill(ZZ.begin(),ZZ.end(),0);    
      std::fill(AA.begin(),AA.end(),0); 
      std::fill(LSPosX.begin(),LSPosX.end(),0.0); 
      std::fill(LSPosY.begin(),LSPosY.end(),0.0); 
      std::fill(LSPosZ.begin(),LSPosZ.end(),0.0); 
      std::fill(LSTheta.begin(),LSTheta.end(),0.0); 
      std::fill(LSPhi.begin(),LSPhi.end(),0.0); 

      std::fill(LSEdep.begin(),LSEdep.end(),0.0);   
      std::fill(Z.begin(),Z.end(),0);    
      std::fill(A.begin(),A.end(),0); 
    }
  }

  TarZ = 0;
  TarA = 0;
  TarEnergy = 0.0;
  TarEdep = 0.0;
  TarPosX = 0.0;
  TarPosY = 0.0; 
  TarPosZ = 0.0;
  TarTheta = 0.0; 
  TarPhi = 0.0;
  TarToF = 0.0;
  Tar_n_e = 0.0;
  Tar_n_theta = 0.0;

  StripEnergy = 0.0;
  StripPosX = 0.0; 
  StripPosY = 0.0; 
  StripPosZ = 0.0; 
  StripTheta = 0.0; 
  StripPhi = 0.0; 

  G4PrimaryVertex* primaryVertex = evt->GetPrimaryVertex();
  G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
  primaryKE = primaryParticle->GetKineticEnergy();

}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event = evt->GetEventID();
  // Target tree
  if(TarEnergy>0.0) {
    // if(TarZ==11) G4cout<<TarEnergy<<" "<<primaryKE<<G4endl;
    analysis->Target(TarZ, TarA, TarEnergy, TarEdep, TarPosX, TarPosY, TarPosZ, TarTheta, TarPhi, StripEnergy, StripPosX, StripPosY, StripPosZ, StripTheta, StripPhi, TarToF, Tar_n_e, Tar_n_theta);
  }
  
  if(LScin)
  {
    for(G4int i = 0; i<LScinNum; i++) 
    {
      // Timing resolution obtained from source spectra with gamma rays
      G4double sigmaT = 1.1;  
      LSToF[i] = G4RandGauss::shoot(LSToF[i], sigmaT);

      // Target - LS coincidence tree
      if(LSEdep[i]>0.01) // || LSEkin[i]>0.0) 
      {  
        if(Z[i]==1 && A[i]==1){
          // G4cout<<event<<"\t"<<i<<"\t"<<primaryKE<<"\t"<<Tar_n_e<<"\t"<<Tar_n_theta<<"\t"<<LSEdep[i]<<G4endl;
          if(LSEdep[i]>0.1) analysis->Target_LS_coinc(TarEnergy, TarTheta, Tar_n_e, Tar_n_theta, LSEkin, LSToF, ZZ, AA, LSPosX, LSPosY, LSPosZ, LSTheta, LSPhi, LSEdep);
        }
        else if(Z[i]==6) {
          if(LSEdep[i]>3.00) analysis->Target_LS_coinc(TarEnergy, TarTheta, Tar_n_e, Tar_n_theta, LSEkin, LSToF, ZZ, AA, LSPosX, LSPosY, LSPosZ, LSTheta, LSPhi, LSEdep);
        }
        else if (Z[i]<1 && A[i]<1){
          if(LSEdep[i]>0.01) analysis->Target_LS_coinc(TarEnergy, TarTheta, Tar_n_e, Tar_n_theta, LSEkin, LSToF, ZZ, AA, LSPosX, LSPosY, LSPosZ, LSTheta, LSPhi, LSEdep);
        } 
      }

      // LS hit tree
      if(LSEkin[i]>0.0) {
        analysis->Reached_LS(LSEkin, LSToF, ZZ, AA, LSPosX, LSPosY, LSPosZ, LSTheta, LSPhi, Tar_n_e, Tar_n_theta);
      }

      // LS detection tree
      if(Z[i]==1 && A[i]==1) //for protons
      {
        if(LSEdep[i]>0.1) {
          analysis->Detected_LS(LSEdep, Z, A);
        }
      }
      else if(Z[i]==6) {  //for carbon 
        if(LSEdep[i]>3.00) {
          analysis->Detected_LS(LSEdep, Z, A);
        }
      }
      //else if (Z[i]<1 && A[i]<1) //for gammas and electrons
      //analysis->Detected_LS(LSEdep, Z, A); 
    }
  }

}   


