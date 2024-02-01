//
// Created by DELL on 2023/12/8.
//

#ifndef G4STARTEventAction_h
#define G4STARTEventAction_h 1
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4THitsMap.hh"
#include "G4STARTSteppingAction.hh"
//class RunAction;
class G4Event;
//class G4THitsMap;
//class G4SDManager;

class G4STARTEventAction : public G4UserEventAction
{
public:
    G4STARTEventAction();
    virtual ~G4STARTEventAction() override ;

    virtual void BeginOfEventAction(const G4Event* anEvent) override ;
    virtual void EndOfEventAction(const G4Event* anEvent) override ;

    //void AddAbs(G4double de, G4double dl);
    //void AddGap(G4double de, G4double dl);

private:
    //RunAction* fRunAction = nullptr;
    //G4double fEdep = 0;
    G4int fEventCounter;
    G4STARTSteppingAction* steppingAction;
   // G4double  fEnergyAbs = 0.;
   // G4double  fEnergyGap = 0.;
   // G4double  fTrackLAbs = 0.;
   // G4double  fTrackLGap = 0.;
    //void WriteToFile(const G4Event* anEvent);
};
/*inline void G4STARTEventAction::AddAbs(G4double de, G4double dl) {
    fEnergyAbs += de;
    fTrackLAbs += dl;
}

inline void G4STARTEventAction::AddGap(G4double de, G4double dl) {
    fEnergyGap += de;
    fTrackLGap += dl;
}*/
#endif //G4START_G4STARTEVENTACTION_HH
