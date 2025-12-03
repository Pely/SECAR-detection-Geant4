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

#ifndef SensitiveDetectorHit_SiMon_h
#define SensitiveDetectorHit_SiMon_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh" // FOR MT

class SensitiveDetectorHit_SiMon : public G4VHit
{
  public:
    SensitiveDetectorHit_SiMon();
    SensitiveDetectorHit_SiMon(const SensitiveDetectorHit_SiMon&);
    virtual ~SensitiveDetectorHit_SiMon();

    // operators
    const SensitiveDetectorHit_SiMon& operator=(const SensitiveDetectorHit_SiMon&);
    G4int operator==(const SensitiveDetectorHit_SiMon&) const;

    inline void *operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    G4int GetID() const { return fId; };

    void SetCounter(G4int count){fcount=count;};
    G4int GetCounter() const {return fcount;};

    void SetEdep(G4double de) { fEdep = de; };
    G4double GetEdep() const  { return fEdep; };

    void SetEkin(G4double de) { fEkin = de; };
    G4double GetEkin() const  { return fEkin; };
    
    // void SetEkin_p(G4double de) { fEkin_p = de; };
    // G4double GetEkin_p() const  { return fEkin_p; };

    void SetBaryonNum(G4int A) { fBaryonNum = A; };
    G4int GetBaryonNumber() const  { return fBaryonNum; };

    void SetParticleCharge(G4int Z) { fCharge = Z; };
    G4int GetParticleCharge() const  { return fCharge; };

    void SetParticleName(G4String name) { fParticleName = name; };
    G4String GetParticleName() const  { return fParticleName; };
     
private:
      G4int      fId, fcount; 
      G4double	 fEdep, fEkin; //, fEkin_p;
      G4double   fCharge, fBaryonNum; 
      G4String   fParticleName;     
};

typedef G4THitsCollection<SensitiveDetectorHit_SiMon> SensitiveDetector_SiMonHitsCollection;

extern G4ThreadLocal G4Allocator<SensitiveDetectorHit_SiMon>* SensitiveDetectorHit_SiMonAllocator;

inline void* SensitiveDetectorHit_SiMon::operator new(size_t)
{
  if(!SensitiveDetectorHit_SiMonAllocator)
      SensitiveDetectorHit_SiMonAllocator = new G4Allocator<SensitiveDetectorHit_SiMon>;
  return (void *) SensitiveDetectorHit_SiMonAllocator->MallocSingle();
}

inline void SensitiveDetectorHit_SiMon::operator delete(void *hit)
{
  SensitiveDetectorHit_SiMonAllocator->FreeSingle((SensitiveDetectorHit_SiMon*) hit);
}

#endif
