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
// Code based on basic example B02


#include "SensitiveDetectorHit_SiMon.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

// THIS IS NECESSARY FOR MT MODE
G4ThreadLocal G4Allocator<SensitiveDetectorHit_SiMon>* SensitiveDetectorHit_SiMonAllocator=0;

SensitiveDetectorHit_SiMon::SensitiveDetectorHit_SiMon()
 : G4VHit(),
   fId(0),
   fEdep(0),
   fEkin(0.)
{}

SensitiveDetectorHit_SiMon::~SensitiveDetectorHit_SiMon() 
{}

SensitiveDetectorHit_SiMon::SensitiveDetectorHit_SiMon(const SensitiveDetectorHit_SiMon& right)
  : G4VHit()
{
  fId         = right.fId;
  fEkin       = right.fEkin;
  fEdep       = right.fEdep;
}

const SensitiveDetectorHit_SiMon& SensitiveDetectorHit_SiMon::operator=(const SensitiveDetectorHit_SiMon& right)
{
  fId         = right.fId;
  fEdep       = right.fEdep;
  fEkin       = right.fEkin;
  return *this;
}

G4int SensitiveDetectorHit_SiMon::operator==(const SensitiveDetectorHit_SiMon& right) const
{
  return ( this == &right ) ? 1 : 0;
}

void SensitiveDetectorHit_SiMon::Draw()
{}

void SensitiveDetectorHit_SiMon::Print()
{
  G4cout<< "HIT: "<< fId << std::setw(6) <<  "Ekin: " <<G4BestUnit(fEkin,"Energy")<< G4endl;
}
//the setw() command sets the space from left to right as to where to write on the terminal
 

