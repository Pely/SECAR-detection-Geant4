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

#include "SteppingAction.hh"
#include "AnalysisManager.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh" 
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"

SteppingAction::SteppingAction(EventAction* eventAction, DetectorConstruction* det, AnalysisManager* pAnalysis)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  detector(det)
{ 
  analysis = pAnalysis; 
}

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4Track* aTrack = aStep -> GetTrack();
  
  G4StepPoint* initialPoint = aStep->GetPreStepPoint();
  G4VPhysicalVolume* volume = initialPoint->GetTouchableHandle()->GetVolume();
  G4String volumeName = volume->GetName();
  G4int Z = aTrack -> GetDynamicParticle() -> GetDefinition() -> GetPDGCharge();
  G4int A = aTrack -> GetDynamicParticle() -> GetDefinition() -> GetBaryonNumber(); 
  // particleName = aTrack->GetDefinition()->GetParticleName();
  // G4StepPoint* finalPoint = aStep->GetPostStepPoint();
  // G4int parentID = aTrack->GetParentID();
  // aTrack->SetTrackStatus(fStopAndKill);
  
  // if(volumeName=="Sphere_phys");
  if(volumeName=="Tar_phys") {
    tarEdep = aStep->GetTotalEnergyDeposit();
    hitTarPos = initialPoint->GetPosition();
    tarDirection = aTrack->GetMomentumDirection();
    // tarEkin = aTrack->GetKineticEnergy();
    tarEkin = aStep->GetPostStepPoint()->GetKineticEnergy();
  }
  if(volumeName=="Strip_phys")
  {
    G4ThreeVector hitStripPos = initialPoint->GetPosition();
    G4ThreeVector stripDirection = aTrack->GetMomentumDirection();
    stripEkin = aStep->GetPostStepPoint()->GetKineticEnergy();
    G4double tarToF = aTrack->GetGlobalTime();
    tarEdep = aStep->GetTotalEnergyDeposit();
    // if(Z==11) G4cout<<tarEkin<<G4endl;
    if(Z==0&&A==1) { // check for neutrons
      theta_n_tar = tarDirection.theta(); ekin_n_tar = tarEkin;
      //ekin_r_tar = 0.0; theta_r_tar = 0.0; 

      //Comment out the following line if doing a reaction and invoke next if statement
      // if(theta_n_tar<0.6)fEventAction->AddTarget(Z, A, tarEkin, hitTarPos, tarDirection, tarToF, ekin_n_tar, theta_n_tar); 
    }
    else {theta_n_tar = 0.0; ekin_n_tar = 0.0; }
    //if(A==58) { // (Z==27 && A==58) store the recoil 
    // ekin_r_tar = tarEkin; theta_r_tar = tarDirection.theta(); 

    //The ekin_n_tar and theta_n_tar data stored here are not useful since they can be distinguished from Z,A 
    //I am keeping this format for now though not efficient
    if(Z!=36) {fEventAction->AddTarget(Z, A, tarEkin, tarEdep, hitTarPos, tarDirection, stripEkin, hitStripPos, stripDirection, tarToF, ekin_n_tar, theta_n_tar); }
  }

  // Flag for turning on/off the Liquid Scintillators. 
  // It is controlled from the header file of the DetectorConstruction 
  if(detector->LScinOn()) 
  { 
    G4int LScinNum = detector->GetLScinNum();
    for(int i=0; i<LScinNum; i++)
    {
      G4double edep=0.0;
      if (volume == detector->GetLScin(i))
      {
        if(initialPoint->GetStepStatus() == fGeomBoundary)
        {
          G4ThreeVector hitPos = initialPoint->GetPosition();
          G4ThreeVector direction = aTrack->GetMomentumDirection();
          G4double ekin = initialPoint->GetKineticEnergy(); 
          G4double tof = aTrack->GetGlobalTime();
          // G4double tof = aTrack->GetLocalTime(); // Use this tof when simulating a source

          if(ekin>0.0001) 
          {
            ekin_LS = ekin; tof_LS = tof; hitPos_LS = hitPos; direction_LS = direction;
            fEventAction->AddLSEntryHit(i, ekin_LS, tof_LS, Z, A, hitPos_LS, direction_LS, ekin_n_tar, theta_n_tar); 
          }
        }

        edep = aStep->GetTotalEnergyDeposit();

        if(edep>0.001) // a very low threshold per step
        {
          fEventAction->AddLSDetected(i,edep,Z,A);
          fEventAction->AddEvent(ekin_r_tar, theta_r_tar, ekin_n_tar, theta_n_tar, i, ekin_LS, tof_LS, Z, A, hitPos_LS, direction_LS, edep);
        }
      }
    }
  }
}

