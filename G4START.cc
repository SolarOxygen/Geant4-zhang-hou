//内置头文件
#include "G4RunManager.hh"
//#include "QGSP_BERT_HP.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4SteppingVerbose.hh"
//#include "G4SystemOfUnits.hh"

//用户编写的头文件
#include "G4STARTDetectorConstruction.hh"
#include "G4STARTPrimaryGeneratorAction.hh"
#include "G4STARTPhysicsList.hh"
#include "G4STARTEventAction.hh"
#include "G4STARTSteppingAction.hh"
//#include "G4STARTRunAction.hh"
int main(int argc,char** argv){
    G4UIExecutive* ui = nullptr;

    if(argc == 1)
    {
        ui = new G4UIExecutive(argc,argv);
    }
    G4int precision = 4;
    G4SteppingVerbose::UseBestUnit(precision);
    //创建G4运行对象
    G4RunManager *runManager = new G4RunManager;

    //G4运行对象参数配置
    //1.建立物理过程
    runManager->SetUserInitialization(new G4STARTDetectorConstruction);
    //2.描述物理过程
    runManager->SetUserInitialization(new G4STARTPhysicsList);
    //3.配置初级事件
    runManager->SetUserAction(new G4STARTPrimaryGeneratorAction);
    runManager->SetUserAction(new G4STARTEventAction);
    //runManager->SetUserAction(new G4STARTRunAction);
    runManager->SetUserAction(new G4STARTSteppingAction);
    //runManager->SetUserAction(new G4STARTRunAction);
    //G4内核初始化
    runManager->Initialize();
    runManager->BeamOn(1);
    auto* visManager = new G4VisExecutive;
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    if ( ! ui ) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }
    else {
        // interactive mode
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    }

    //启动一次


    //释放内存
    delete visManager;
    delete runManager;

    return 0;
}