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

#include "PrimaryGeneratorAction.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(AnalysisManager* pAnalysis)
{
    analysis = pAnalysis;  

    if(source)
    {
        fParticleGun  = new G4ParticleGun(1);
        fParticleGun->SetParticleEnergy(0*eV);
        fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-31.*mm));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.)); 
    }
    if(pn_neutrons)
    {
        n=1;
        fParticleGun  = new G4ParticleGun(1);
        fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
        fParticleGun->SetParticleDefinition(particle);
    }
    if(pn_recoils)
    {
        n=1;
        fParticleGun  = new G4ParticleGun(1);
        fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-31.*mm));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
        fParticleGun->SetParticleEnergy(0*MeV);
    }
    else gps = new G4GeneralParticleSource();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    if(source) delete fParticleGun; 
    else delete gps;
}	

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if(source)
    {
        fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-30.*mm));

        // Distribution uniform in solid angle
        G4double cosTheta = 2*G4UniformRand() - 1.; 
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        G4double phi = twopi*G4UniformRand(); 
        G4double vx = sinTheta*std::cos(phi),
                 vy = sinTheta*std::sin(phi),
                 vz = cosTheta;

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else if(pn_neutrons )
    {  

        G4double energy, theta;
        if(n==1) OpenExternalFile();
        
        input_file >> theta >> energy;

        G4double theta_n = 0.0;
        G4double energy_n = 0.0;
        if(energy>0) energy_n = CLHEP::RandFlat::shoot(energy-0.002, energy+0.002);
        else energy_n = CLHEP::RandFlat::shoot(energy, energy+0.002);
        if(theta>0) theta_n = CLHEP::RandFlat::shoot(theta-0.002, theta+0.002);
        else theta_n = theta+0.000001;


        fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
        // Uniform spherical distribution
        // Converting the spherical coordinates (polar and azimuthal angles) to Cartesian coordinates. 
        G4double phi = CLHEP::RandFlat::shoot(0.0,360.);
        G4double vx = sin(theta_n*deg)*cos(phi*deg);
        G4double vy = sin(theta_n*deg)*sin(phi*deg);
        G4double vz = cos(theta_n*deg);

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
        fParticleGun->SetParticleEnergy(energy_n*MeV);  
        fParticleGun->GeneratePrimaryVertex(anEvent);
        n++;
        //if(n==nmax) {n=0; input_file.close();} //This is the line where the file ends
        if(input_file.eof()){ input_file.close(); n=1; } 
    }
    else if(pn_recoils)
    {  
        if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {
            G4int Z = 27 , A = 58;
            G4double ionCharge = 0.*eplus;
            G4double excitEnergy = 0.*keV;
            G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
            fParticleGun->SetParticleCharge(ionCharge);
            fParticleGun->SetParticleDefinition(ion);
        }

        G4double energy, theta;
        if(n==1) OpenExternalFile();
        
        input_file >> theta >> energy;

        G4double theta_r = 0.0;
        G4double energy_r = 0.0;
        if(energy>0) energy_r = CLHEP::RandFlat::shoot(energy-0.002, energy+0.002);
        else energy_r = CLHEP::RandFlat::shoot(energy, energy+0.002);
        if(theta>0) theta_r = CLHEP::RandFlat::shoot(theta-0.002, theta+0.002);
        else theta_r = theta+0.000001;


        fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-20*cm));
        // Uniform spherical distribution
        // Converting the spherical coordinates (polar and azimuthal angles) to Cartesian coordinates. 
        G4double phi = CLHEP::RandFlat::shoot(0.0,360.);
        G4double vx = sin(theta_r*deg)*cos(phi*deg);
        G4double vy = sin(theta_r*deg)*sin(phi*deg);
        G4double vz = cos(theta_r*deg);

        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
        fParticleGun->SetParticleEnergy(energy_r*MeV);  
        fParticleGun->GeneratePrimaryVertex(anEvent);
        n++;
        //if(n==nmax) {n=0; input_file.close();} //This is the line where the file ends
        if(input_file.eof()){ input_file.close(); n=1; } 
    }
    else gps->GeneratePrimaryVertex(anEvent);

    
}

void PrimaryGeneratorAction::OpenExternalFile()
{   
    // Open neutron distribution file
    if(pn_neutrons)
    input_file.open("../neutronInput/pn_distributions/geant_input/input_n_0_219.txt",std::ios::in);
    if(pn_recoils)
    input_file.open("input/input_rec_3_211.txt",std::ios::in);
    
    // !!!!! Important !!!!! Make sure the file doesn't have an empty line at the end
    // If there is remove it, it will create warnings while running !!!!

    // if (input_file.is_open() && n==1) G4cout << "Neutrons are genereted using external file!! "<< G4endl;
    // else G4cout << "File not found!! "<< G4endl;
    
}
