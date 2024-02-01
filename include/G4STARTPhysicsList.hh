//
// Created by DELL on 2023/12/7.
//

#ifndef G4STARTPhysicsList_h
#define G4STARTPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
class G4STARTPhysicsList:public G4VModularPhysicsList
{
public:
    G4STARTPhysicsList();
    ~G4STARTPhysicsList() override ;

    void SetCuts() override ;
};
#endif //G4START_G4STARTPHYSICSLISTS_HH
