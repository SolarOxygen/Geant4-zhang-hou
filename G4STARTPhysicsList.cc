//
// Created by DELL on 2023/12/7.
//
#include "G4STARTPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scintillation.hh"
#include "G4ScintillationTrackInformation.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
G4STARTPhysicsList::G4STARTPhysicsList():G4VModularPhysicsList()
{
    SetVerboseLevel(1);
    defaultCutValue = 0.1*mm;

    RegisterPhysics(new G4DecayPhysics());
    //RegisterPhysics(new G4DecayPhysics());

    // EM physics
    RegisterPhysics(new G4EmStandardPhysics());

    //G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    //RegisterPhysics(opticalPhysics);

    RegisterPhysics(new G4RadioactiveDecayPhysics());
    //RegisterPhysics(new G4EmLivermorePolarizedPhysics());
}
G4STARTPhysicsList::~G4STARTPhysicsList() {}

void G4STARTPhysicsList::SetCuts()
{
    G4VUserPhysicsList::SetCuts();
}