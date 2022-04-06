#ifndef Scintillator_H
#define Scintillator_H 1
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4String.hh"

class Scintillator
{
    public:
    Scintillator();
    void Place(G4RotationMatrix *pRot, 
                        G4ThreeVector &tlate, 
                        const G4String &pName, 
                        G4LogicalVolume *pMotherLogical,  
                        G4int pCopyNo);
    void ConstructSDandField();
    
    private:
    G4Material* MakePlactic();
    G4Material* plasticScint;
    void ConstructFront();
    void ConstructRear();
    void ConstructSingleSD(G4String name, G4LogicalVolume* logVol);
    G4LogicalVolume* scintFrontLog;
    G4LogicalVolume* scintRear2hLog;
    G4LogicalVolume* scintRear4hLog;
    G4LogicalVolume* scintRear6hLog;
    G4LogicalVolume* scintRear8hLog;
    G4LogicalVolume* scintRear10hLog;
    G4LogicalVolume* scintRear12hLog;
};

#endif
