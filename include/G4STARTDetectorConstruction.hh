#ifndef G4STARTDetectorConstruction_h
#define G4STARTDetectorConstruction_h 1


#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "tls.hh"

class G4Material;
class G4MaterialPropertiesTable;
class G4STARTElectricFieldSetup;

class G4STARTDetectorConstruction : public G4VUserDetectorConstruction{
public:
    G4STARTDetectorConstruction();
    ~G4STARTDetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void DefineMaterials();
    virtual void ConstructSDandField() override;
    //void SetMaxStep (G4double );
protected:
    G4Material* diamond;
    //G4MaterialPropertiesTable* diamondMPT;
    //G4Material* fPb;
    G4Material* fW;
    G4Material* fSi;
    G4Material* fAl;

private:
    G4Cache<G4STARTElectricFieldSetup*> fEmFieldSetup;
};
#endif