#ifndef G4STARTPrimaryGeneratorAction_h
#define G4STARTPrimaryGeneratorAction_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"



class G4Event;

class G4STARTPrimaryGeneratorAction :public G4VUserPrimaryGeneratorAction{

public:
    //构造函数
    G4STARTPrimaryGeneratorAction();
    //析构函数
    ~G4STARTPrimaryGeneratorAction() override;

    // 事件生成的函数
    virtual void GeneratePrimaries(G4Event* anEvent) override;
private:
    //声明粒子枪
    G4ParticleGun* particleGun;
    G4double x;
    G4double y;
};
#endif