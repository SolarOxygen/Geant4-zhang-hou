//内置头文件
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4UserLimits.hh"
#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4AutoDelete.hh"
#include "G4UniformMagField.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
//#include "G4DormandPrince745.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


#include "G4SDManager.hh"
#include "G4PSEnergyDeposit.hh"
//#include "G4PSNofSecondary.hh"
//#include "G4PSNofStep.hh"

//用户编写的头文件
#include "G4STARTDetectorConstruction.hh"
#include "G4STARTElectricFieldSetup.hh"
G4STARTDetectorConstruction::G4STARTDetectorConstruction() {
    DefineMaterials();
}
G4STARTDetectorConstruction::~G4STARTDetectorConstruction() {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    delete sdManager;
}

void G4STARTDetectorConstruction::DefineMaterials() {
    /*G4double a;
    G4double z;
    G4double density;
    G4Element* fC= new G4Element("C","C", z = 6.,a = 12.0107* g/mole );
    diamond =new G4Material("diamond",density =3.50 * g/cm3,1);
    diamond->AddElement(fC,1);
    G4MaterialPropertiesTable* diamondMPT = new G4MaterialPropertiesTable();

    G4double photonEnergy[3] = {2.0 * eV, 3.0 * eV, 4.0 * eV}; // Example energy levels
    G4double refractiveIndex[3] = {2.42, 2.40, 2.38};
    G4double absorptionLength[3] = {100.0 * mm, 50.0 * mm, 30.0 * mm};
    G4double scatteringLength[] = {0.05, 0.07, 0.1};
    diamondMPT->AddProperty("RINDEX",photonEnergy,refractiveIndex,3);
    diamondMPT->AddProperty("ABSLENGTH",photonEnergy,absorptionLength,3);
    diamondMPT->AddProperty("RAYLEIGH", photonEnergy, scatteringLength, 3);

    diamond->SetMaterialPropertiesTable(diamondMPT);
    delete diamondMPT;

    fPb = new G4Material("Pb",82.,207.2*g/mole,11.3437*g/cm3);
    G4MaterialPropertiesTable* PbMPT = new G4MaterialPropertiesTable();

// Example data (replace with actual values):
    G4double photonEnergy2[] = {1.0, 2.0, 3.0};  // Replace with actual energy values
    G4double refractiveIndex2[] = {2.5, 2.3, 2.1};  // Replace with actual refractive index values
    G4double absorptionLength2[] = {10.0, 15.0, 20.0};  // Replace with actual absorption length values
    G4double scatteringLength2[] = {5.0, 7.0, 10.0};  // Replace with actual scattering length values
    G4double efficiency[] = {0.8, 0.9, 1.0};  // Replace with actual efficiency values

    PbMPT->AddProperty("RINDEX", photonEnergy2, refractiveIndex2, 3);
    PbMPT->AddProperty("ABSLENGTH", photonEnergy2, absorptionLength2, 3);
    PbMPT->AddProperty("RAYLEIGH", photonEnergy2, scatteringLength2, 3);
    PbMPT->AddProperty("EFFICIENCY", photonEnergy2, efficiency, 3);

    fPb->SetMaterialPropertiesTable(PbMPT);
    delete PbMPT;*/
    G4NistManager* nist = G4NistManager::Instance();
    G4String materialName = "Diamond";
    G4double density = 3.52 * g/cm3;  // Density of diamond
    G4double z = 6.0;  // Atomic number of carbon
    G4double a = 12.01 * g/mole;  // Atomic mass of carbon
    diamond = new G4Material(materialName, density, 1, kStateSolid);

    diamond->AddElement(nist->FindOrBuildElement("C"), 1);

    fW = nist->FindOrBuildMaterial("G4_W");
    fAl = nist->FindOrBuildMaterial("G4_Al");
    fSi = nist->FindOrBuildMaterial("G4_Si");

    //delete nist;
}
G4VPhysicalVolume* G4STARTDetectorConstruction::Construct() {
    //创建材料查找表
    G4NistManager* nist2 = G4NistManager::Instance();
    //创建world
    //创建world立方体的几何体
    G4Material* vacuumMaterial = nist2->FindOrBuildMaterial("G4_Galactic");

    // If not available, create a custom vacuum material
    if (!vacuumMaterial) {
        G4String name = "Custom_Vacuum";
        G4double density = 1.0e-25 * g/cm3;  // Very low density
        G4double temperature = 0.1 * kelvin; // A low temperature
        G4double pressure = 1.0e-19 * pascal; // A low pressure

        vacuumMaterial = new G4Material(name, density, 0, kStateGas, temperature, pressure);
    }
    G4Box* solid_world = new G4Box("world",                 //实体名称
                                    5*cm,                  //x轴长度
                                    5*cm,                  //y轴
                                    550*mm                   //z轴
                                   );

    //创建world逻辑体
    G4LogicalVolume* logic_world = new G4LogicalVolume(solid_world,                             //待填充的几何体
                                                        vacuumMaterial,    //填充材料
                                                         "world"                                //逻辑体名称
                                                         );
    //创建world物理实体
    G4PVPlacement* phy_world = new G4PVPlacement(0,                         //旋转
                                                 G4ThreeVector(0,0,0),      //坐标位置
                                                 logic_world,               //待摆放的逻辑体
                                                 "world",                   //物理实体名称
                                                 0,                         //母体
                                                 false,                     //是否布尔运算
                                                 0                          //编号
                                                 );

    //创建G4的靶材料
    //创建靶材料几何体

    G4Tubs* tub_one = new G4Tubs("Tub_one",
                                 10*mm,
                                 12*mm,
                                 42.5*mm,
                                 0.*deg,
                                 360.*deg);
    //创建靶材料逻辑体
    G4LogicalVolume* logic_tub_one = new G4LogicalVolume(tub_one,                             //待填充的几何体
                                                       fSi,    //填充材料
                                                       "Tub_one"                                //逻辑体名称
                                                        );

    //创建靶材料物理实体
    G4PVPlacement* phy_tub_one = new G4PVPlacement(0,                         //旋转
                                                 G4ThreeVector(0,0,-507.5*mm),      //坐标位置
                                                  logic_tub_one,               //待摆放的逻辑体
                                                 "Tub_one",                   //物理实体名称
                                                 logic_world,                         //母体
                                                 false,                     //是否布尔运算
                                                 1                          //编号
                                                 );
           //shape 3

    G4Tubs* tub_two = new G4Tubs("Tub_two",
                                 10*mm,
                                 12*mm,
                                 47.5*mm,
                                 0.*deg,
                                 360.*deg
                                 );

    //创建W靶材料逻辑体
    G4LogicalVolume* logic_tub_two = new G4LogicalVolume(tub_two,                             //待填充的几何体
                                                        fAl,    //填充材料
                                                        "Tub_two"                                //逻辑体名称
    );
    //创建靶材料物理实体
    G4PVPlacement* phy_tub_two = new G4PVPlacement(0,                         //旋转
                                                  G4ThreeVector(0,0,-417.5*mm),      //坐标位置
                                                  logic_tub_two,               //待摆放的逻辑体
                                                  "Tub_two",                   //物理实体名称
                                                  logic_world,                         //母体
                                                  false,                     //是否布尔运算
                                                  1                          //编号
    );
    G4Cons* cons_one = new G4Cons("Cons_one",
                                  10*mm,
                                  12*mm,
                                  17.5*mm,
                                  19.5*mm,
                                  237.5*mm,
                                  0.*deg,
                                  360.*deg
                                  );
    G4LogicalVolume* logic_cons_one = new G4LogicalVolume(cons_one,                             //待填充的几何体
                                                         fAl,    //填充材料
                                                         "Cons_one"                                //逻辑体名称
    );
    G4PVPlacement* phy_cons_one = new G4PVPlacement(0,                         //旋转
                                                   G4ThreeVector(0,0,-132.5*mm),      //坐标位置
                                                   logic_cons_one,               //待摆放的逻辑体
                                                   "Cons_one",                   //物理实体名称
                                                   logic_world,                         //母体
                                                   false,                     //是否布尔运算
                                                   1                          //编号
    );
    G4Cons* cons_two = new G4Cons("Cons_two",
                                  0.1*mm,
                                  0.2*mm,
                                  10*mm,
                                  12*mm,
                                  237.5*mm,
                                  0.*deg,
                                  360.*deg
    );
    G4LogicalVolume* logic_cons_two = new G4LogicalVolume(cons_two,                             //待填充的几何体
                                                          fAl,    //填充材料
                                                          "Cons_two"                                //逻辑体名称
    );
    G4PVPlacement* phy_cons_two = new G4PVPlacement(0,                         //旋转
                                                    G4ThreeVector(0,0,-132.5*mm),      //坐标位置
                                                    logic_cons_two,               //待摆放的逻辑体
                                                    "Cons_two",                   //物理实体名称
                                                    logic_world,                         //母体
                                                    false,                     //是否布尔运算
                                                    1                          //编号
    );
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

    // declare crystal as a MultiFunctionalDetector scorer
    //
    G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
    G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
    G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
    cryst->RegisterPrimitive(primitiv1);
    logic_tub_one->SetSensitiveDetector(cryst);

    auto visAttributes = new G4VisAttributes(G4Colour(0,1.0,0));
    logic_tub_two->SetVisAttributes(visAttributes);
    auto visAttributes2 = new G4VisAttributes(G4Colour(1.0,0,0));
    logic_tub_one->SetVisAttributes(visAttributes2);
    auto visAttributes3 = new G4VisAttributes(G4Colour(0,0,1.0));
    logic_cons_one->SetVisAttributes(visAttributes3);
    auto visAttributes4 = new G4VisAttributes(G4Colour(1.0,1.0,0));
    logic_cons_two->SetVisAttributes(visAttributes4);
    /*G4UniformElectricField(G4ThreeVector(0., 0., 0.));
    G4ElectricField* pEMfield;
    pEMfield = new G4UniformElectricField(
            G4ThreeVector(0.,0.,200000*volt/m));
    G4FieldManager* fieldManager= G4TransportationManager::GetTransportationManager()->
            GetFieldManager();
// Set this field to the global field manager
    fieldManager->SetDetectorField( pEMfield );
    G4MagneticField *magField;
    magField = new G4UniformMagField(G4ThreeVector(0.,3.0*kilogauss,0.));
    fieldManager->SetDetectorField( magField );
    fieldManager->CreateChordFinder( magField );*/

    return phy_world;
}
void G4STARTDetectorConstruction::ConstructSDandField() {

    //delete fieldManager;
    /*if (!fEmFieldSetup.Get()) {
        G4STARTElectricFieldSetup* fieldSetup = new G4STARTElectricFieldSetup();
        G4AutoDelete::Register(fieldSetup); //Kernel will delete the messenger
        fEmFieldSetup.Put(fieldSetup);
    }*/
}
