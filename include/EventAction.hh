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

#ifndef EventAction_h
#define EventAction_h 1

// General Geant4 libraries
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"


//C++ Stuff
#include <vector>
#include <fstream>
#include <iostream>

class AnalysisManager;
class DetectorConstruction;

class EventAction : public G4UserEventAction
{
public:

  EventAction(AnalysisManager*, DetectorConstruction*);
  ~EventAction();

  virtual void  BeginOfEventAction(const G4Event* evt);
  virtual void    EndOfEventAction(const G4Event* evt);


  
  void AddTarget(G4int z, G4int a, G4double e, G4double de, G4ThreeVector pos, G4ThreeVector ang, G4double e_str, G4ThreeVector pos_str, G4ThreeVector ang_str, G4double dt, G4double e_n, G4double th_n)
  {
    TarZ = z; 
    TarA = a; 
    TarEnergy = e;
    TarEdep += de;
    TarPosX = pos.x(); 
    TarPosY = pos.y(); 
    TarPosZ = pos.z(); 
    TarTheta = ang.theta(); 
    TarPhi = ang.phi(); 
    TarToF = dt;
    Tar_n_e = e_n;
    Tar_n_theta = th_n;

    StripEnergy = e_str;
    StripPosX = pos_str.x(); 
    StripPosY = pos_str.y(); 
    StripPosZ = pos_str.z(); 
    StripTheta = ang_str.theta(); 
    StripPhi = ang_str.phi(); 

  };

  void AddLSDetected(G4int i, G4double de, G4int z, G4int a)
  { 
    LSEdep[i] += de; 
    Z[i] = z; 
    A[i] = a; 
  };

  void AddLSEntryHit(G4int i, G4double e, G4double dt, G4int z, G4int a, G4ThreeVector pos, G4ThreeVector ang, G4double e_nt, G4double th_nt)
  { 
    LSEkin[i] = e; 
    LSToF[i] = dt; 
    ZZ[i] = z; 
    AA[i] = a; 
    LSPosX[i] = pos.x(); 
    LSPosY[i] = pos.y(); 
    LSPosZ[i] = pos.z(); 
    LSTheta[i] = ang.theta(); 
    LSPhi[i] = ang.phi();
    Tar_n_e = e_nt;
    Tar_n_theta = th_nt;
  };

  void AddEvent(G4double e_rt, G4double th_rt, G4double e_nt, G4double th_nt, G4int i, G4double e, G4double t, G4int z, G4int a, G4ThreeVector pos, G4ThreeVector ang, G4double de)
  {
    TarEnergy = e_rt;
    TarTheta = th_rt; 
    Tar_n_e = e_nt;
    Tar_n_theta = th_nt;

    LSEkin[i] = e; 
    LSToF[i] = t; 
    ZZ[i] = z; 
    AA[i] = a; 
    LSPosX[i] = pos.x(); 
    LSPosY[i] = pos.y(); 
    LSPosZ[i] = pos.z(); 
    LSTheta[i] = ang.theta(); 
    LSPhi[i] = ang.phi();

    LSEdep[i] += de;
  };

private:
      
  AnalysisManager* analysis;
  DetectorConstruction* detector;

  G4double TarPosX, TarPosY, TarPosZ, TarTheta, TarPhi, TarEnergy, TarEdep, TarToF, Tar_n_e, Tar_n_theta;
  G4int TarZ, TarA; 
  G4double StripPosX, StripPosY, StripPosZ, StripTheta, StripPhi, StripEnergy;

  G4double primaryKE;

  typedef std::vector<G4double> DVec;
  DVec LSPosX, LSPosY, LSPosZ, LSTheta, LSPhi; 
  DVec LSEkin, LSEdep, LSToF;
  
  typedef std::vector<G4int> IVec;
  IVec Z, A, ZZ, AA;

  G4int LScinNum;
  G4bool LScin;
};
#endif

    
