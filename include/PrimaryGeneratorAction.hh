
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "TF1.h"
#include <vector>

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
     PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

    void GeneratePrimaries(G4Event* anEvent);
    enum class TransitionType { FF, GT };

  private:
    G4ParticleGun* particleGun;	
	void SetUp( void );
	void SetUpBetaAssymetry();
	void SetUpGammaAssymetry();
	//direction
	void GenerateAsymDirection(G4ThreeVector* direction, TF1* asymDistr);		
	TF1* FindAsymDistrFunc(double asymFactor, double velocity, double polarisation);
	double FindAsymmertyFactor(TransitionType transitionType, int deltaI=0, int initialSpin = 0);
    double FindVelocity(double energyinMeV);
    G4ThreeVector GenerateIsotropicDirection( G4double thetaMin = 0,
                                              G4double thetaMax = M_PI,
                                              G4double phiMin = 0,
                                              G4double phiMax = 2.*M_PI);
                                              
                                              
    TF1* betaAsymDistr;
    TF1* gammaAsymDistr;
    double betaEnergy;
     
};

#endif // PrimaryGeneratorAction_h
