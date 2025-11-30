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

#ifndef SensitiveDetectorHit_h
#define SensitiveDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh" // FOR MT

class SensitiveDetectorHit : public G4VHit
{
  public:
    SensitiveDetectorHit();
    SensitiveDetectorHit(const SensitiveDetectorHit&);
    virtual ~SensitiveDetectorHit();

    // operators
    const SensitiveDetectorHit& operator=(const SensitiveDetectorHit&);
    G4int operator==(const SensitiveDetectorHit&) const;

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

    void SetPosition(CLHEP::Hep3Vector pos) { fPos_x = pos.x(); fPos_y = pos.y(); fPos_z = pos.z(); };
    G4double GetPositionX() const  { return fPos_x; };
    G4double GetPositionY() const  { return fPos_y; };
    G4double GetPositionZ() const  { return fPos_z; };

    void SetMomentum(CLHEP::Hep3Vector mom) { fMomentum_x = mom.x(); fMomentum_y = mom.y(); fMomentum_z = mom.z(); };
    G4double GetMomentumX() const  { return fMomentum_x; };
    G4double GetMomentumY() const  { return fMomentum_y; };
    G4double GetMomentumZ() const  { return fMomentum_z; };

    void SetAngle(CLHEP::Hep3Vector angle) { fTheta = angle.theta(); fPhi = angle.phi(); };
    G4double GetAngleTheta() const  { return fTheta; };
    G4double GetAnglePhi() const    { return fPhi; };
    
    void SetBaryonNum(G4int A) { fBaryonNum = A; };
    G4int GetBaryonNumber() const  { return fBaryonNum; };

    void SetParticleCharge(G4int Z) { fCharge = Z; };
    G4int GetParticleCharge() const  { return fCharge; };

    void SetParticleName(G4String name) { fParticleName = name; };
    G4String GetParticleName() const  { return fParticleName; };
     
private:
      G4int      fId, fcount; 
      G4double	 fEdep, fEkin;
      G4double   fPos_x, fPos_y, fPos_z;
      G4double   fMomentum_x, fMomentum_y, fMomentum_z;
      G4double   fTheta, fPhi;
      G4double   fCharge, fBaryonNum; 
      G4String   fParticleName;     
};

typedef G4THitsCollection<SensitiveDetectorHit> SensitiveDetectorHitsCollection;

extern G4ThreadLocal G4Allocator<SensitiveDetectorHit>* SensitiveDetectorHitAllocator;

inline void* SensitiveDetectorHit::operator new(size_t)
{
  if(!SensitiveDetectorHitAllocator)
      SensitiveDetectorHitAllocator = new G4Allocator<SensitiveDetectorHit>;
  return (void *) SensitiveDetectorHitAllocator->MallocSingle();
}

inline void SensitiveDetectorHit::operator delete(void *hit)
{
  SensitiveDetectorHitAllocator->FreeSingle((SensitiveDetectorHit*) hit);
}

#endif
