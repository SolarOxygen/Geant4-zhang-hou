#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4STARTPrimaryGeneratorAction.hh"
#include "G4UniformRandPool.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "Randomize.hh"


//构造函数
G4STARTPrimaryGeneratorAction::G4STARTPrimaryGeneratorAction(){
    //定义一个粒子源，粒子枪一次打出一个粒子
    particleGun = new G4ParticleGun(1);

    //生成粒子源查找表
    G4ParticleTable* table = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle
            = table->FindParticle(particleName="e-");

    //设置打出粒子的类型，计划打出中子
    particleGun->SetParticleDefinition(particle);

    //设置粒子的能量
    G4double energy=0;
    const std::vector<G4double> energyValues = {100.0*MeV,150.0*MeV,200.0*MeV,250.0*MeV,500.0*MeV,1000.0*MeV};
    const std::vector<G4double> energyWeights ={0.1,0.1,0.1,0.2,0.2,0.3};
    G4double randNum = G4UniformRand();
    G4double cumulativeWeight = 0.0;
    for (size_t i = 0; i < energyValues.size(); ++i) {
        cumulativeWeight += energyWeights[i];
        if (randNum <= cumulativeWeight) {
            energy=energyValues[i];
            break;
        }
    }
    particleGun->SetParticleEnergy(1000.0*MeV);
    //particleGun1->SetParticleEnergy(1*eV);

    //设置粒子源位置
    /*G4double radius = 0.1*m;
    G4double theta = (G4UniformRand()-0.5) * 2.0 * pi;
    G4double x = radius * std::cos(theta);
    G4double y = radius * std::sin(theta);
    particleGun->SetParticlePosition(G4ThreeVector(x,y,-45*cm));
    //particleGun->SetParticlePosition(G4ThreeVector(0,0,-45*cm));*/
    //设置粒子打出方向;
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));


}
//析构函数
G4STARTPrimaryGeneratorAction::~G4STARTPrimaryGeneratorAction(){
    delete particleGun;

}

//事件生成函数
void G4STARTPrimaryGeneratorAction:: GeneratePrimaries(G4Event* e){

    G4double radius = 0.5*cm;
    //G4double theta = (G4UniformRand()-0.5) * 2.0 * pi;
    x = radius * (G4UniformRand()-0.5)*2;
    y = radius * (G4UniformRand()-0.5)*2;
    std::ofstream outFile("particle_position_diamond.txt", std::ios::app);
    outFile << x/cm <<" "<< y/cm << "\n";
    outFile.close();
    particleGun->SetParticlePosition(G4ThreeVector(x,y,-0.3*cm));
    //发射粒子
    particleGun->GeneratePrimaryVertex(e);
}