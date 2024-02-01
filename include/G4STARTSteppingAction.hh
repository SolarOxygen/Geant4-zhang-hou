//
// Created by DELL on 2023/12/9.
//

#ifndef G4STARTSteppingAction_h
#define G4STARTSteppingAction_h 1
#include "G4UserSteppingAction.hh"
#include <map>
#include "G4String.hh"
//class EventAction;
class G4STARTSteppingAction:public G4UserSteppingAction
{
public:
    G4STARTSteppingAction();
    ~G4STARTSteppingAction() override;
    void UserSteppingAction(const G4Step* step) override;
    const std::map<G4String, G4int>& GetSecondaryParticleCounts() const;
    void ResetCounters() ;
    G4int GetElectronHolePairs() const;
    void ResetElectronHolePairs();

    //G4double GetCrossSection() const;
private:
    //G4double totalComptonEvents;
    //EventAction* fEventAction = nullptr;
    std::map<G4String, G4int> secondaryParticleCounts;
    G4int electronHolePairs;
};

#endif //G4START_G4STARTSTEPPINGACTION_HH
