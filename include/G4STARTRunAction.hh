//
// Created by DELL on 2023/12/8.
//

#ifndef G4STARTRunAction_h
#define G4STARTRunAction_h 1
#include "G4UserRunAction.hh"
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include <fstream>
#include "G4STARTSteppingAction.hh"
class G4STARTRunAction : public G4UserRunAction {
public:
    G4STARTRunAction();
    virtual ~G4STARTRunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);

private:
    //std::ofstream outputFile;  // Output file stream
    G4STARTSteppingAction* steppingAction;
};

#endif //G4START_G4STARTRUNACTION_HH
