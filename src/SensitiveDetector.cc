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
// Author: Pelagia Tsintari, tsint1p@cmich.edu 
//
// Code based on the advanced example radioprotection

#include "SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4INCLRandom.hh"

SensitiveDetector::SensitiveDetector(const G4String& name, const G4String& hitsCollectionName, AnalysisManager* analysis_manager) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
  analysis = analysis_manager;
  count = 0;
}

SensitiveDetector::~SensitiveDetector() 
{}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new SensitiveDetectorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{  
  G4Track* aTrack = aStep->GetTrack();

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ThreeVector momentum = aTrack->GetMomentumDirection();
  G4ThreeVector position = aTrack->GetPosition();
  G4int charge = aTrack->GetDefinition()->GetPDGCharge();  
  G4int baryon = aTrack->GetDefinition()->GetBaryonNumber();
  G4String particleName = aTrack->GetDefinition()->GetParticleName();


  if (ekin==0.) return false;
  G4String volumeName = aStep -> GetTrack() -> GetVolume() -> GetName();
  if(volumeName != "Detector_phys") return false;  

  SensitiveDetectorHit* newHit = new SensitiveDetectorHit();

  newHit -> SetEdep(edep);
  newHit -> SetEkin(ekin);
  newHit -> SetPosition(position);
  newHit -> SetAngle(momentum);
  newHit -> SetParticleName(particleName);
  newHit -> SetParticleCharge(charge);
  newHit -> SetBaryonNum(baryon);
   
  //if(charge==18&&baryon==40)
  // if(charge==26&&baryon==58)
  // {
  //   count++;
  //   if(count==10001)
  //     {newHit->SetCounter(0);}
  //   else
  //     {newHit->SetCounter(count);}
  // }

  fHitsCollection -> insert( newHit );

  return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{

  //Initialisation of total energy deposition per event to zero
  G4double totalEdep = 0.0*MeV;
  G4double Ekin = 0.0*MeV;
  
  G4int Z = 0;
  G4int A = 0;
  G4double pos_x = 0;
  G4double pos_y = 0;
  G4double theta = 0;
  G4double phi   = 0;
  G4String particleName = "";
  
  // G4bool neutron;
  // G4bool recoil;
  // G4double Ekin_R = 0.0*MeV;
  // G4double Ekin_n = 0.0*MeV;
  // G4double theta_n = 0;
  // G4double posn_x = 0;
  // G4double posn_y = 0;
  // G4int coincidence=0;
  // G4int beam = 0;


  G4int NbHits = fHitsCollection->entries();

  for (G4int i=0;i<NbHits;i++)
  {
    G4double edep = (*fHitsCollection)[i]->GetEdep();
    totalEdep = totalEdep + edep;
    particleName  = (*fHitsCollection)[0]->GetParticleName();
    Ekin  = (*fHitsCollection)[0]->GetEkin();
    A     = (*fHitsCollection)[0]->GetBaryonNumber();
    Z     = (*fHitsCollection)[0]->GetParticleCharge();
    pos_x = (*fHitsCollection)[0]->GetPositionX();
    pos_y = (*fHitsCollection)[0]->GetPositionY();
    phi   = (*fHitsCollection)[0]->GetAnglePhi(); 
    count = (*fHitsCollection)[i]->GetCounter();

    //if(Z==18&A==40)
    // if(Z==26&A==58)
    // { 
    //   if(count==10000) beam=0; 
    //   else beam=1;  //Flag for not storing the beam per hit
    // }

    //if(Z==19&A==40)
    // if(Z==27&A==58)
    // {
    //   recoil=true;
    //   Ekin_R=Ekin;
    // }

    // if(particleName=="neutron")
    // {
    //   neutron=true;
    //   Ekin_n=Ekin;
    //   posn_x=pos_x;
    //   posn_y=pos_y;
    //   theta_n=57.296*atan(sqrt((posn_x*posn_x)+(posn_y*posn_y))/500);
    // }
    
    // if((neutron==true)&(recoil==true))
    // {
    //   coincidence=1;
    //   neutron=false;
    //   recoil=false;
    // }
    // else coincidence=0;

  }
  //if(beam==0)
  G4double sigmaE = 0.2; //σ = 150 keV -> FWHM = 150*2.355
  totalEdep= G4RandGauss::shoot(totalEdep, sigmaE);   

  theta=57.296*atan(sqrt((pos_x*pos_x)+(pos_y*pos_y))/700);
  // if(Z!=26)
  if(totalEdep>0)
  {
    if(Z==11) analysis->Detector_DSSD(totalEdep, Ekin, A, Z, theta, phi*57.296, pos_x, pos_y);
  }//, coincidence, Ekin_R, Ekin_n, posn_x, posn_y, theta_n); 

}





