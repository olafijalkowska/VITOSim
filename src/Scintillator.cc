
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh" 
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4MultiFunctionalDetector.hh"
#include "Scintillator.hh"
#include "CADMesh.hh"
#include <iostream>

Scintillator::Scintillator()
{
    plasticScint = MakePlactic();
	ConstructFront();
	ConstructRear();
}


G4Material* Scintillator::MakePlactic()
{
    G4NistManager* man=G4NistManager::Instance();
    G4Element* H = man->FindOrBuildElement("H");
    G4Element* C = man->FindOrBuildElement("C");
    G4Material* EJ200 = new G4Material("EJ200", 1.023*g/cm3, 2);
    EJ200->AddElement(H, 517);
    EJ200->AddElement(C, 469);
    return EJ200;

}


void Scintillator::ConstructFront()
{
    CADMesh* frontMesh = new CADMesh("../stl/scint-FrontDet.stl",
                                    1*mm, G4ThreeVector(0*cm, 0*cm, 0*cm), 0);
                                    
    G4VSolid* frontSolid =  frontMesh->TessellatedMesh(); 
    scintFrontLog = new G4LogicalVolume(frontSolid, plasticScint, "scintFrontLog");
    G4VisAttributes* scintFrontVisAtt = new G4VisAttributes(G4Colour::Yellow());
	scintFrontVisAtt->SetForceAuxEdgeVisible(true);// Can see outline when drawn with lines
	scintFrontVisAtt->SetForceSolid(true);
	scintFrontLog->SetVisAttributes(scintFrontVisAtt);
}

void Scintillator::ConstructRear()
{
    CADMesh* rearMesh2h = new CADMesh("../stl/Rear-2h.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);
    /*CADMesh* rearMesh2h = new CADMesh("../stl/scintil.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);    */                            
                                    
                                                                       
    G4VSolid* rearSoli2h =  rearMesh2h->TessellatedMesh(); 
    scintRear2hLog = new G4LogicalVolume(rearSoli2h, plasticScint, "rearSoli2h");
    
    CADMesh* rearMesh4h = new CADMesh("../stl/Rear-4h.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);                                    
    G4VSolid* rearSoli4h =  rearMesh4h->TessellatedMesh(); 
    scintRear4hLog = new G4LogicalVolume(rearSoli4h, plasticScint, "rearSoli4h");
    
    CADMesh* rearMesh6h = new CADMesh("../stl/Rear-6h.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);                                    
    G4VSolid* rearSoli6h =  rearMesh6h->TessellatedMesh(); 
    scintRear6hLog = new G4LogicalVolume(rearSoli6h, plasticScint, "rearSoli6h");
    
    
    CADMesh* rearMesh8h = new CADMesh("../stl/Rear-8h.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);                                    
    G4VSolid* rearSoli8h =  rearMesh8h->TessellatedMesh(); 
    scintRear8hLog = new G4LogicalVolume(rearSoli8h, plasticScint, "rearSoli8h");
    
    CADMesh* rearMesh10h = new CADMesh("../stl/Rear-10.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);                                    
    G4VSolid* rearSoli10h =  rearMesh10h->TessellatedMesh(); 
    scintRear10hLog = new G4LogicalVolume(rearSoli10h, plasticScint, "rearSoli10h");
    
    CADMesh* rearMesh12h = new CADMesh("../stl/Rear-12.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);                                    
    G4VSolid* rearSoli12h =  rearMesh12h->TessellatedMesh(); 
    scintRear12hLog = new G4LogicalVolume(rearSoli12h, plasticScint, "rearSoli12h");        
            
    G4VisAttributes* scintRearVisAtt = new G4VisAttributes(G4Colour::Green());
	scintRearVisAtt->SetForceAuxEdgeVisible(true);// Can see outline when drawn with lines
	scintRearVisAtt->SetForceSolid(true);
	scintRear2hLog->SetVisAttributes(scintRearVisAtt);
	scintRear4hLog->SetVisAttributes(scintRearVisAtt);
	scintRear6hLog->SetVisAttributes(scintRearVisAtt);
	scintRear8hLog->SetVisAttributes(scintRearVisAtt);
	scintRear10hLog->SetVisAttributes(scintRearVisAtt);
	scintRear12hLog->SetVisAttributes(scintRearVisAtt);
	
	/*G4VisAttributes* scintFrontVisAtt = new G4VisAttributes(G4Colour::Yellow());
	scintFrontVisAtt->SetForceAuxEdgeVisible(true);// Can see outline when drawn with lines
	scintFrontVisAtt->SetForceSolid(true);
	scintRear4hLog->SetVisAttributes(scintFrontVisAtt);
	scintRear6hLog->SetVisAttributes(scintFrontVisAtt);
	scintRear8hLog->SetVisAttributes(scintFrontVisAtt);
	scintRear10hLog->SetVisAttributes(scintFrontVisAtt);
	scintRear12hLog->SetVisAttributes(scintFrontVisAtt);*/
	
}




void Scintillator::Place(G4RotationMatrix *pRot, 
                        G4ThreeVector &tlate, 
                        const G4String &pName, 
                        G4LogicalVolume *pMotherLogical,  
                        G4int pCopyNo)
{
    new G4PVPlacement(pRot, tlate, scintFrontLog, pName+"Front",  pMotherLogical, 0,pCopyNo);
    new G4PVPlacement(pRot, tlate, scintRear2hLog, pName+"Rear2h",  pMotherLogical, 0,pCopyNo);
    new G4PVPlacement(pRot, tlate, scintRear4hLog, pName+"Rear4h",  pMotherLogical, 0,pCopyNo);
    new G4PVPlacement(pRot, tlate, scintRear6hLog, pName+"Rear6h",  pMotherLogical, 0,pCopyNo);
    new G4PVPlacement(pRot, tlate, scintRear8hLog, pName+"Rear8h",  pMotherLogical, 0,pCopyNo);
    new G4PVPlacement(pRot, tlate, scintRear10hLog, pName+"Rear10h",  pMotherLogical, 0,pCopyNo);
    new G4PVPlacement(pRot, tlate, scintRear12hLog, pName+"Rear12h",  pMotherLogical, 0,pCopyNo);
}

void Scintillator::ConstructSingleSD(G4String name, G4LogicalVolume* logicVol)
{
    G4MultiFunctionalDetector* detector = new G4MultiFunctionalDetector(name);
    G4int depth = 0;
    G4VPrimitiveScorer* energyDepScorer = new G4PSEnergyDeposit("eDep",depth);
    detector->RegisterPrimitive(energyDepScorer);
    logicVol->SetSensitiveDetector(detector);
    G4SDManager* SDmanager = G4SDManager::GetSDMpointer();
    SDmanager->AddNewDetector(detector);

}

void Scintillator::ConstructSDandField()
{
	ConstructSingleSD("scintFrontSD", scintFrontLog);
	/*G4MultiFunctionalDetector* detector = new G4MultiFunctionalDetector("scintFrontSD");
    G4int depth = 0;
    G4VPrimitiveScorer* energyDepScorer = new G4PSEnergyDeposit("eDep",depth);
    detector->RegisterPrimitive(energyDepScorer);
    scintFrontLog->SetSensitiveDetector(detector);
    G4SDManager* SDmanager = G4SDManager::GetSDMpointer();
    SDmanager->AddNewDetector(detector);*/
	
	ConstructSingleSD("scintRear2hSD", scintRear2hLog);
	ConstructSingleSD("scintRear2hSD", scintRear2hLog);
	ConstructSingleSD("scintRear4hSD", scintRear4hLog);
	ConstructSingleSD("scintRear6hSD", scintRear6hLog);
	ConstructSingleSD("scintRear8hSD", scintRear8hLog);
	ConstructSingleSD("scintRear10hSD", scintRear10hLog);
	ConstructSingleSD("scintRear12hSD", scintRear12hLog);
}
