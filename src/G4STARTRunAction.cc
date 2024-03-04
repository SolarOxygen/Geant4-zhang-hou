//
// Created by DELL on 2023/12/8.
//
#include "G4STARTRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4RootAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4STARTSteppingAction.hh"

G4STARTRunAction::G4STARTRunAction() : G4UserRunAction(),steppingAction(){}

G4STARTRunAction::~G4STARTRunAction() {
    // Close the output file stream

}

void G4STARTRunAction::BeginOfRunAction(const G4Run* run) {
    // Open the output file for writing
    //outputFile.open("particle_in_W.txt");
    steppingAction = dynamic_cast<G4STARTSteppingAction*>(G4EventManager::GetEventManager()->GetUserSteppingAction());
}

void G4STARTRunAction::EndOfRunAction(const G4Run* run) {
    // Close the output file
    if (steppingAction) {
        // Write final state particle counts to the file
        const std::map<G4String, G4int>& counts = steppingAction->GetSecondaryParticleCounts();
        //std::ofstream outFile("particle_in_W.txt", std::ios::app);
        for (const auto& pair : counts) {
            std::string name = pair.first;
            char filename[20];
            std::sprintf(filename,"particle_in_Si_%s.txt",name);
            std::ofstream outFile(filename, std::ios::app);
            outFile <<pair.second << "\n";
            outFile.close();
            //G4cout << pair.first<<" " << pair.second<<"\n";
        }

        steppingAction->ResetCounters();
    }
}