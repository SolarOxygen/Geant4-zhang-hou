//
// Created by DELL on 2023/12/8.
//
#include "G4STARTEventAction.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofSecondary.hh"
#include "G4AnalysisManager.hh"
#include <fstream>
#include "G4STARTSteppingAction.hh"
#include "G4EventManager.hh"
//#include "G4RootAnalysisManager.hh"
G4STARTEventAction::G4STARTEventAction() : fEventCounter(0),steppingAction(){}

G4STARTEventAction::~G4STARTEventAction() {}

void G4STARTEventAction::BeginOfEventAction(const G4Event *anEvent) {
    steppingAction = dynamic_cast<G4STARTSteppingAction*>(G4EventManager::GetEventManager()->GetUserSteppingAction());
}

void G4STARTEventAction::EndOfEventAction(const G4Event *anEvent) {
    G4SDManager *sdManager = G4SDManager::GetSDMpointer();
    G4HCofThisEvent* HCE = anEvent->GetHCofThisEvent();
    if(!HCE)
        return;;
    G4int edepCollectionID = sdManager->GetCollectionID("edep");
    G4THitsMap<G4double> *edepHitsMap =(G4THitsMap<G4double>*)(HCE->GetHC(edepCollectionID));
    std::map<G4int,G4double*>::iterator itr;

    //G4int nSecCollectionID = sdManager->GetCollectionID("NumberOfSecondaries");
    //G4THitsMap<G4int> *nSecHitsMap = dynamic_cast<G4THitsMap<G4int> *>(anEvent->GetHCofThisEvent()->GetHC(
          //  nSecCollectionID));
    if (edepHitsMap ) {
        // Open a text file for writing
        std::ofstream outFile("particle_info_diamond.txt", std::ios::app);
        // Write information for each hit
        int i=0;
        itr = edepHitsMap->GetMap()->begin();
        if(itr==edepHitsMap->GetMap()->end()&&fEventCounter!=0)
        {
            outFile << 0 <<" "<< fEventCounter <<"\n";
        }
        for (itr = edepHitsMap->GetMap()->begin();itr!=edepHitsMap->GetMap()->end();itr++) {
            if(fEventCounter==0){
                i++;
                break;
            }
            G4double edep = *(itr->second);
            //int *numSecondaries = (*nSecHitsMap)[i];
            outFile << edep/keV <<" "<< fEventCounter <<"\n";
            //outFile << "Number of Secondaries: " << numSecondaries << "\n";
            i++;
        }

        // Close the text file
        outFile.close();

        // Increment event counter
        fEventCounter++;
    }
    if (steppingAction) {
        // Write final state particle counts to the file
        const std::map<G4String, G4int>& counts = steppingAction->GetSecondaryParticleCounts();
        //std::ofstream outFile("particle_in_W.txt", std::ios::app);
        int k=0;
        int pa[3]={0,0,0};
        for (const auto& pair : counts) {
            std::string name = pair.first;
            char filename[20];
            std::sprintf(filename,"particle_in_W_%s.txt",name);
            std::ofstream outFile(filename, std::ios::app);
            if(name=="W182"||name=="W183"||name=="W184"||name=="W186")
            {

            }
            else
            {
                if(name=="gamma")
                {
                    pa[0]++;
                }
                else if(name=="e-")
                {
                    pa[1]++;
                }
                else if(name=="e+")
                {
                    pa[2]++;
                }
                k++;
            }
            outFile <<pair.second <<" "<<fEventCounter-1 << "\n";
            outFile.close();

            //k++;
            //G4cout << pair.first<<" " << pair.second<<"\n";
        }
        if (fEventCounter!=1 && k<3) {
            if(pa[0]==0) {
                std::ofstream outFile("particle_in_W_gamma.txt", std::ios::app);
                outFile << 0 << " " << fEventCounter - 1 << "\n";
                outFile.close();
            }
            if(pa[1]==0) {
                std::ofstream outFile1("particle_in_W_e-.txt", std::ios::app);
                outFile1 << 0 << " " << fEventCounter - 1 << "\n";
                outFile1.close();
            }
            if(pa[2]==0) {
                std::ofstream outFile2("particle_in_W_e+.txt", std::ios::app);
                outFile2 << 0 << " " << fEventCounter - 1 << "\n";
                outFile2.close();
            }
        }
        steppingAction->ResetCounters();
        G4int electronHolePairs = steppingAction->GetElectronHolePairs();
        std::ofstream outFile3("electron_pair_in_diamond.txt", std::ios::app);
        if(fEventCounter-1 !=0) {
            outFile3 << electronHolePairs << " " << fEventCounter - 1 << "\n";
        }
        outFile3.close();
        steppingAction->ResetElectronHolePairs();
    }
}