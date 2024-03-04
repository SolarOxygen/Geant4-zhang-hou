//
// Created by DELL on 2023/12/9.
//
#include "G4STARTSteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include <map>
#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
G4STARTSteppingAction::G4STARTSteppingAction() :G4UserSteppingAction(){

}

G4STARTSteppingAction::~G4STARTSteppingAction() {

}

void G4STARTSteppingAction::UserSteppingAction(const G4Step* step)
{
    const G4Track* track = step->GetTrack();
    const G4ParticleDefinition* particle = track->GetParticleDefinition();
        //if(particle == G4Gamma::Definition()){
        const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
        if(process){
            G4ThreeVector  momentumDirection = track->GetMomentumDirection();
            G4double energy = track->GetKineticEnergy();
            if(energy == 0*MeV) {}
            else{
                G4String particleName = particle->GetParticleName();
                std::ofstream outFile("particle_info2.txt", std::ios::app);
                outFile << "\n Energy of " <<particleName << "="<< energy / MeV << "MeV";
                outFile << "\n Momentum Direction: " << momentumDirection << "\n";
                outFile << "-----------------\n";

            }
        }
   // }

    if(track->GetParentID() > 0 ){
        G4String particleName = particle->GetParticleName();
        //secondaryParticleCounts[particleName]=0;
        secondaryParticleCounts[particleName]++;
    }
    if (track) {
        // Check if the particle is in the diamond detector
        if (track->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() == "G4_Si") {
            // Check if an ionization event occurred (e.g., an electron-hole pair is created)
            if (step->GetTotalEnergyDeposit() > 5.5*eV) {
                // Increment the count for each electron-hole pair
                int num = step->GetTotalEnergyDeposit()/(13.1*eV);
                electronHolePairs+= num;
            }
        }
    }
}
const std::map<G4String, G4int>& G4STARTSteppingAction::GetSecondaryParticleCounts() const {
    return secondaryParticleCounts;
}
void G4STARTSteppingAction:: ResetCounters(){
    secondaryParticleCounts.clear();
}
G4int G4STARTSteppingAction::GetElectronHolePairs() const {
    return electronHolePairs;
}
void G4STARTSteppingAction::ResetElectronHolePairs() {
    electronHolePairs = 0;
}