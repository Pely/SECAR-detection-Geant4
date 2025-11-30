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
// Code based on basic example B02


#include "SensitiveDetectorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

// THIS IS NECESSARY FOR MT MODE
G4ThreadLocal G4Allocator<SensitiveDetectorHit>* SensitiveDetectorHitAllocator=0;

SensitiveDetectorHit::SensitiveDetectorHit()
 : G4VHit(),
   fId(0),
   fEdep(0),
   fEkin(0.),
   fPos_x(0.),
   fPos_y(0.),
   fPos_z(0.),
   fMomentum_x(0.),
   fMomentum_y(0.),
   fMomentum_z(0.),
   fTheta(0.),
   fPhi(0.)
{}

SensitiveDetectorHit::~SensitiveDetectorHit() 
{}

SensitiveDetectorHit::SensitiveDetectorHit(const SensitiveDetectorHit& right)
  : G4VHit()
{
  fId         = right.fId;
  fEkin       = right.fEkin;
  fEdep       = right.fEdep;
  fPos_x      = right.fPos_x;
  fPos_y      = right.fPos_y;
  fPos_z      = right.fPos_z;
  fMomentum_x = right.fMomentum_x;
  fMomentum_y = right.fMomentum_y;
  fMomentum_z = right.fMomentum_z;
  fTheta      = right.fTheta;
  fPhi        = right.fPhi;
}

const SensitiveDetectorHit& SensitiveDetectorHit::operator=(const SensitiveDetectorHit& right)
{
  fId         = right.fId;
  fEdep       = right.fEdep;
  fEkin       = right.fEkin;
  fPos_x      = right.fPos_x;
  fPos_y      = right.fPos_y;
  fPos_z      = right.fPos_z;
  fMomentum_x = right.fMomentum_x;
  fMomentum_y = right.fMomentum_y;
  fMomentum_z = right.fMomentum_z;
  fTheta      = right.fTheta;
  fPhi        = right.fPhi;
  return *this;
}

G4int SensitiveDetectorHit::operator==(const SensitiveDetectorHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

void SensitiveDetectorHit::Draw()
{}

void SensitiveDetectorHit::Print()
{
  G4cout<< "HIT: "<< fId << std::setw(6) <<  "Ekin: " <<G4BestUnit(fEkin,"Energy")<< G4endl;
}
//the setw() command sets the space from left to right as to where to write on the terminal
 

