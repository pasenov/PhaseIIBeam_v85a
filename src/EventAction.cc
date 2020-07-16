//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: EventAction.cc 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"
#include <iomanip>

#include <stdio.h>
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* DA, RunAction* RA)
:G4UserEventAction(), fDetectorconstruction(DA), fRunAction(RA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
 //geometry parameters
 posZAl=fDetectorconstruction->GetPosZAl();
 posZScint1=fDetectorconstruction->GetPosZScint1();
 posZScint2=fDetectorconstruction->GetPosZScint2();
 posZBPIX12=fDetectorconstruction->GetPosZBPIX12();
 posZBPIX34=fDetectorconstruction->GetPosZBPIX34();
 posZ1=fDetectorconstruction->GetPosZ1();
 posZ2=fDetectorconstruction->GetPosZ2();
 posZDUT=fDetectorconstruction->GetPosZDUT();
 Det1SizeZ=fDetectorconstruction->GetSize1();
 Det2SizeZ=fDetectorconstruction->GetSize2();
 Dist=fDetectorconstruction->GetDist();
 Strip1Depth=fDetectorconstruction->GetStrip1Depth();
 Strip1Length=fDetectorconstruction->GetStrip1Length();
 Strip2Depth=fDetectorconstruction->GetStrip2Depth();
 Strip2Length=fDetectorconstruction->GetStrip2Length();
 StripDist=fDetectorconstruction->GetStripDist();
 StripWidth=fDetectorconstruction->GetStripWidth();
 StripPitch = StripDist + StripWidth;
 posEndArm1=-(fDetectorconstruction->Getpos_EndArm1Abs());
 posBeginningArm2=fDetectorconstruction->Getpos_BeginningArm2Abs();
 BPIXSizeZ=fDetectorconstruction->GetSizeBPIX();
 pixelX=fDetectorconstruction->GetPixelPitchX();
 pixelY=fDetectorconstruction->GetPixelPitchY();
 pixelDepth=fDetectorconstruction->GetPixelDepth();
 XangleDUT=fDetectorconstruction->GetDUTangleX();
 XangleBPIX=fDetectorconstruction->GetBPIXangleX();
 YangleBPIX=fDetectorconstruction->GetBPIXangleY();
 ElField1=fDetectorconstruction->GetElField1();
 ElField2=fDetectorconstruction->GetElField2();
 totalNbStrips = 1016;
 ROChor = 52*pixelX + 2*pixelX;
 ROCvert = 80*pixelY + pixelY;
 fPixThreshold = 10000;
 fStripThreshold = 5000;
 fPixSigma = 1400;
 fStripSigma = 1500;
 iRandLim = 50;
 stubLim = 100*um;

 //G4cout 
     //<< "\n DetectorConstruction parameters assigned. "
     //<< G4endl;

 // initialisation per event
 fEnergyDeposit1  =  fEnergyDeposit2 = fEnergyDepositL1  =  fEnergyDepositL2 = fEnergyDepositROCL1  =  fEnergyDepositROCL2 = fEnergyDepositSc1 = fEnergyDepositSc2 = fEnergyDepositAir = fEnergyDepositAl = fEnergySecondary1 = fEnergySecondary2 = fEnergyTertiary1 = fEnergyTertiary2 = 0.;
 fEnergyEntSensor1 = fEnergyExSensor1 = fEnergyEntSensor2 = fEnergyExSensor2 = fEnergyEntPixelLayer1 = fEnergyExPixelLayer1 = fEnergyEntPixelLayer2 = fEnergyExPixelLayer2 = fEnergyEntScintillator1 = fEnergyExScintillator1 = fEnergyEntScintillator2 = fEnergyExScintillator2 = fEnergyExAl = 0.;
 fTimeEntranceScint1 = fTimeEntranceScint2 = 0.;
 fNb2SElectrons = 0;
 for (G4int j = 0; j <= 1016; j = j + 1) {
    fEnergyStrip1a[j] = 0.;
    fEnergyStrip1b[j] = 0.;
    fEnergyStrip2a[j] = 0.;
    fEnergyStrip2b[j] = 0.;
    fWeightStrip1a[j] = 0.;
    fWeightStrip1b[j] = 0.;
    fWeightStrip2a[j] = 0.;
    fWeightStrip2b[j] = 0.;
    fNbHitsStrip1a[j] = 0;
    fNbHitsStrip1b[j] = 0;
    fNbHitsStrip2a[j] = 0;
    fNbHitsStrip2b[j] = 0;
    fChargeStrip1a[j] = 0;
    fChargeStrip1b[j] = 0;
    fChargeStrip2a[j] = 0;
    fChargeStrip2b[j] = 0;
    fChargeStripMod1a[j] = 0;
    fChargeStripMod1b[j] = 0;
    fChargeStripMod2a[j] = 0;
    fChargeStripMod2b[j] = 0;
    fHitStrip1a[j] = 0;
    fHitStrip1b[j] = 0;
    fHitStrip2a[j] = 0;
    fHitStrip2b[j] = 0;
 }

 totalNbHitStrips = 0;

 for (G4int je1 = 0; je1 <= 2; je1 = je1 + 1) {
    fHitsMultiplicityPix[je1] = 0;
    for (G4int je2 = 0; je2 <= 16; je2 = je2 + 1) {
 	for (G4int je3 = 0; je3 <= 80; je3 = je3 + 1) {
	   for (G4int je4 = 0; je4 <= 52; je4 = je4 + 1) {
		fEnergyPixel[je1][je2][je3][je4] = 0.;
		fChargePixel[je1][je2][je3][je4] = 0;
		fChargePixelMod[je1][je2][je3][je4] = 0;
		fNbHitsPixel[je1][je2][je3][je4] = 0;
 	   }
 	}
    }
 }

 for (G4int je5 = 0; je5 <= 16; je5 = je5 + 1) {
    for (G4int je6 = 0; je6 <= 80; je6 = je6 + 1) {
       for (G4int je7 = 0; je7 <= 52; je7 = je7 + 1) {
          fWeightPixel1[je5][je6][je7] = 0.;
          fWeightPixel2[je5][je6][je7] = 0.;
       }
    }
 }

 fChargeStrip1 = fChargeStrip2 = 0;
 fChargePixel1 = fChargePixel2 = 0;
 fCS1a = fCS1b = fCS2a = fCS2b = 0;
 fChargePix1 = fChargePix2 = 0;
 fChargeQuad = 0;
 fChargePixMod = 0;

 fPrimaryTrackLength = fSecondaryTrackLength = fSecondaryDetTrackLength = fTertiaryTrackLength = 0.;
 fSecondaryxPolarization = fSecondaryyPolarization = fSecondaryzPolarization = fTertiaryxPolarization = fTertiaryyPolarization = fTertiaryzPolarization = 0.;
 fHitsMultiplicity1b = fHitsMultiplicity2b = 0;
 fHitSensor1 = fHitSensor2 = 0;
 fHitPixelDet1 = fHitPixelDet2 = 0;

 fHitMatchFoundInThisEvent = fStubFoundInThisEvent = 0;

 fMomDir1x = fMomDir1y = fMomDir1z =  fMomDir2x = fMomDir2y = fMomDir2z = 0.;

 mod = row = col = -10;

 for (G4int k = 0; k < 5; k = k + 1)  {
    for (G4int i = 0; i < 417; i = i + 1)  {
       for (G4int j = 0; j < 161; j = j + 1)  {
	  fClusterOccupancy[k][i][j] = 0;
       }
    }
 }

 fPointPix1ent = fPointPix1entReal = fPointPix1mid = fPointPix1ex = fPointDet1ent = fPointDet1entReal = fPointDet1mid = fPointDet1ex = fPointDet2ent = fPointDet2entReal = fPointDet2mid = fPointDet2ex = fPointPix2ent = fPointPix2entReal = fPointPix2mid = fPointPix2ex = fPoint20Scint = fPointExitAl = fPointExitScint1 = fPointExitScint2 = fPointExitPixel12 = fPointExitPixel34 = fPointExitSensor1 =  fPointExitSensor2 = G4ThreeVector(0, 0, 0);

 fTrack1 = fTrack2 = 0.;

 dB = dC = dBAD = dCAD = dAD = 0.;
 B1Bx = B1By = C1Cx = C1Cy = 0.;

 fDiff1x = fDiff1y = fDiff2x = fDiff2y = 0.;

 xBA = xBD = xCA = xCD = xAD = xBAD = xCAD = B1 = C1 = G4ThreeVector(0, 0, 0);

 //G4cout 
     //<< "\n End of BeginningOfEventAction. "
     //<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 // get event ID
 G4int evtNb = evt->GetEventID();
 //G4cout 
    //<< "\n Event ID = " 
    //<< evtNb
    //<< G4endl;

 fRunAction->AddEnergyDeposit1(fEnergyDeposit1);
 fRunAction->AddEnergyDeposit2(fEnergyDeposit2);

 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 analysisManager->SetVerboseLevel(1);
 analysisManager->FillH1(1, fEnergyDeposit1);
 analysisManager->FillH1(2, fEnergyDeposit2);
 analysisManager->FillH1(3, fEnergySecondary1);
 analysisManager->FillH1(4, fEnergySecondary2);
 analysisManager->FillH1(5, fEnergyTertiary1);
 analysisManager->FillH1(6, fEnergyTertiary2);
 analysisManager->FillH1(7, fEnergyDeposit1+fEnergySecondary1);
 analysisManager->FillH1(8, fEnergyDeposit2+fEnergySecondary2);
 analysisManager->FillH1(21, fSecondaryxPolarization);
 analysisManager->FillH1(22, fSecondaryyPolarization);
 analysisManager->FillH1(23, fSecondaryzPolarization);
 analysisManager->FillH1(24, fTertiaryxPolarization);
 analysisManager->FillH1(25, fTertiaryyPolarization);
 analysisManager->FillH1(26, fTertiaryzPolarization);
 analysisManager->FillH1(27, fPrimaryTrackLength);
 analysisManager->FillH1(28, fSecondaryTrackLength);
 analysisManager->FillH1(29, fTertiaryTrackLength);
 analysisManager->FillH1(30, fSecondaryDetTrackLength);
 analysisManager->FillH1(31, fTrack2);
 analysisManager->FillH1(212, fEnergyDepositSc1);
 analysisManager->FillH1(213, fEnergyDepositSc2);
 analysisManager->FillH1(214, fEnergyDepositAir);
 analysisManager->FillH1(215, fEnergyDepositAl);
 analysisManager->FillH1(216, fEnergyDeposit1);
 analysisManager->FillH1(217, fEnergyDeposit2);
 analysisManager->FillH1(218, fEnergyDepositL1);
 analysisManager->FillH1(219, fEnergyDepositL2);
 analysisManager->FillH1(220, fEnergyExAl);
 analysisManager->FillH1(221, fEnergyEntScintillator1);
 analysisManager->FillH1(222, fEnergyExScintillator1);
 analysisManager->FillH1(223, fEnergyEntPixelLayer1);
 analysisManager->FillH1(224, fEnergyExPixelLayer1);
 analysisManager->FillH1(225, fEnergyEntSensor1);
 analysisManager->FillH1(226, fEnergyExSensor1);
 analysisManager->FillH1(227, fEnergyEntSensor2);
 analysisManager->FillH1(228, fEnergyExSensor2);
 analysisManager->FillH1(229, fEnergyEntPixelLayer2);
 analysisManager->FillH1(230, fEnergyExPixelLayer2);
 analysisManager->FillH1(231, fEnergyEntScintillator2);
 analysisManager->FillH1(232, fEnergyExScintillator2);
 analysisManager->FillH1(233, fEnergyDepositROCL1);
 analysisManager->FillH1(234, fEnergyDepositROCL2);
 analysisManager->FillH1(235, fTimeEntranceScint1);
 analysisManager->FillH1(236, fTimeEntranceScint2);
 analysisManager->FillH1(237, fTimeEntranceScint2-fTimeEntranceScint1);
 analysisManager->FillH1(238, fNb2SElectrons);

 /*G4cout
    << "\n fNb2SElectrons filled"
    << fNb2SElectrons
    << G4endl;*/

 for (G4int i = 52; i < 60; i = i + 1) {
   analysisManager->FillH1(i, fEnergyStrip1a[i+453]);
 }
 for (G4int p = 60; p < 68; p = p + 1) {
   analysisManager->FillH1(p, fEnergyStrip1b[p+445]);
 }
 for (G4int l = 68; l < 76; l = l + 1) {
   analysisManager->FillH1(l, fEnergyStrip2a[l+437]);
 }
 for (G4int n = 76; n < 84; n = n + 1) {
   analysisManager->FillH1(n, fEnergyStrip2b[n+429]);
 }

 // fill ntuple and charge per strip histograms
 for (G4int k = 0; k < 1017; k++) {
   analysisManager->FillNtupleDColumn(k, fEnergyStrip1a[k]);
   analysisManager->FillNtupleDColumn(k+1017, fEnergyStrip1b[k]);
   analysisManager->FillNtupleDColumn(k+2*1017, fEnergyStrip2a[k]);
   analysisManager->FillNtupleDColumn(k+3*1017, fEnergyStrip2b[k]);

   //charge collected per strip in number of electrons [e]
   fChargeStrip1a[k] = (fEnergyStrip1a[k])/(3.67*eV);
   fChargeStrip1b[k] = (fEnergyStrip1b[k])/(3.67*eV);
   fChargeStrip2a[k] = (fEnergyStrip2a[k])/(3.67*eV);
   fChargeStrip2b[k] = (fEnergyStrip2b[k])/(3.67*eV);

   fChargeStripMod1a[k] = (fEnergyStrip1a[k])/(3.67*eV);
   fChargeStripMod1b[k] = (fEnergyStrip1b[k])/(3.67*eV);
   fChargeStripMod2a[k] = (fEnergyStrip2a[k])/(3.67*eV);
   fChargeStripMod2b[k] = (fEnergyStrip2b[k])/(3.67*eV);

   fCS1a = fChargeStrip1a[k];
   //G4cout 
     //<< "\n fCS1a = "
     //<< fCS1a;
   fCS1b = fChargeStrip1b[k];
   //G4cout 
     //<< "\n fCS1b = "
     //<< fCS1b;
   fCS2a = fChargeStrip2a[k];
   //G4cout 
     //<< "\n fCS2a = "
     //<< fCS2a;
   fCS2b = fChargeStrip2b[k];
   //G4cout 
     //<< "\n fCS2b = "
     //<< fCS2b;
 
   if (fCS1a >= fStripThreshold)  {
      fRunAction->AddNbHitsStrip1a(k);
      fHitSensor1 = 1;
      fHitStrip1a[k] = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 1 Row: a Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip1a advanced"
        //<< G4endl;

      fChargeStrip1 = fChargeStrip1 + fStripThreshold;
      totalNbHitStrips++;
   }
   if (fCS1b >= fStripThreshold)  {
      fRunAction->AddNbHitsStrip1b(k);
      fHitSensor1 = 1;
      //fHitsMultiplicity1b ++;
      fHitStrip1b[k] = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 1 Row: b Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip1b advanced"
        //<< G4endl;

      fChargeStrip1 = fChargeStrip1 + fStripThreshold;
      totalNbHitStrips++;
   }
   if (fCS2a >= fStripThreshold)  {
      fRunAction->AddNbHitsStrip2a(k);
      fHitSensor2 = 1;
      fHitStrip2a[k] = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 2 Row: a Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip2a advanced"
        //<< G4endl;

      fChargeStrip2 = fChargeStrip2 + fStripThreshold;
      totalNbHitStrips++;
   }
   if (fCS2b >= fStripThreshold)  {
      fRunAction->AddNbHitsStrip2b(k);
      fHitSensor2 = 1;
      //fHitsMultiplicity2b ++;
      fHitStrip2b[k] = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 2 Row: b Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip2b advanced"
        //<< G4endl;

      fChargeStrip2 = fChargeStrip2 + fStripThreshold;
      totalNbHitStrips++;
   }
 }

 // charge sharing for strips
 for (G4int k = 0; k < 1017; k++) {
   fRandStrip1 = G4RandGauss::shoot((fChargeStrip1a[k])/2,fStripSigma);
   fRandStrip2 = G4RandGauss::shoot((fChargeStrip1b[k])/2,fStripSigma);
   fRandStrip3 = G4RandGauss::shoot((fChargeStrip2a[k])/2,fStripSigma);
   fRandStrip4 = G4RandGauss::shoot((fChargeStrip2b[k])/2,fStripSigma);
   fRandStrip5 = G4RandGauss::shoot((fChargeStrip1a[k])/2,fStripSigma);
   fRandStrip6 = G4RandGauss::shoot((fChargeStrip1b[k])/2,fStripSigma);
   fRandStrip7 = G4RandGauss::shoot((fChargeStrip2a[k])/2,fStripSigma);
   fRandStrip8 = G4RandGauss::shoot((fChargeStrip2b[k])/2,fStripSigma);

   iSecretStrip1 = rand() % 99;
   iSecretStrip2 = rand() % 99;
   iSecretStrip3 = rand() % 99;
   iSecretStrip4 = rand() % 99;
   iSecretStrip5 = rand() % 99;
   iSecretStrip6 = rand() % 99;
   iSecretStrip7 = rand() % 99;
   iSecretStrip8 = rand() % 99;

   if (k > 1)  {
      if (iSecretStrip1 < iRandLim)  {
        fChargeStripMod1a[k-1]+= (fChargeStrip1a[k])/2 + fRandStrip1;
      }
      if (iSecretStrip2 < iRandLim)  {
        fChargeStripMod1b[k-1]+= (fChargeStrip1b[k])/2 + fRandStrip2;
      }
      if (iSecretStrip3 < iRandLim)  {
        fChargeStripMod2a[k-1]+= (fChargeStrip2a[k])/2 + fRandStrip3;
      }
      if (iSecretStrip4 < iRandLim)  {
        fChargeStripMod2b[k-1]+= (fChargeStrip2b[k])/2 + fRandStrip4;
      }
   }

   if (k <1016)  {
      if (iSecretStrip5 < iRandLim)  {
        fChargeStripMod1a[k+1]+= (fChargeStrip1a[k])/2 + fRandStrip5;
      }
      if (iSecretStrip6 < iRandLim)  {
        fChargeStripMod1b[k+1]+= (fChargeStrip1b[k])/2 + fRandStrip6;
      }
      if (iSecretStrip7 < iRandLim)  {
        fChargeStripMod2a[k+1]+= (fChargeStrip2a[k])/2 + fRandStrip7;
      }
      if (iSecretStrip8 < iRandLim)  {
        fChargeStripMod2b[k+1]+= (fChargeStrip2b[k])/2 + fRandStrip8;
      }
   }
 }
 for (G4int k = 0; k < 1017; k++) {
   fCS1a = fChargeStripMod1a[k];
   //G4cout 
     //<< "\n fCS1a = "
     //<< fCS1a;
   fCS1b = fChargeStripMod1b[k];
   /*G4cout 
     << "\n k = "
     << k
     << "\n fCS1b = "
     << fCS1b
     << G4endl;*/
   fCS2a = fChargeStripMod2a[k];
   //G4cout 
     //<< "\n fCS2a = "
     //<< fCS2a;
   fCS2b = fChargeStripMod2b[k];
   /*G4cout 
     << "\n k = "
     << k
     << "\n fCS2b = "
     << fCS2b
     << G4endl;*/
 
   if (fCS1a >= fStripThreshold)  {
   }
   if (fCS1b >= fStripThreshold)  {
      fHitsMultiplicity1b ++;
   }
   if (fCS2a >= fStripThreshold)  {
   }
   if (fCS2b >= fStripThreshold)  {
      fHitsMultiplicity2b ++;
   }
 }


 // check for hit matches in 2S sensors
 for (G4int k = 0; k < 1017; k++) {
   if (fHitStrip1a[k] == 1)  {
     if (k == 1)  {
       if ((fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1) || (fHitStrip2a[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 2)  {
       if ((fHitStrip2a[k-1] == 1) ||(fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1) || (fHitStrip2a[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if ((k >= 3) && (k <= 1014))  {
       if ((fHitStrip2a[k-2] == 1) || (fHitStrip2a[k-1] == 1) || (fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1) || (fHitStrip2a[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1015)  {
       if ((fHitStrip2a[k-2] == 1) || (fHitStrip2a[k-1] == 1) || (fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1016)  {
       if ((fHitStrip2a[k-2] == 1) || (fHitStrip2a[k-1] == 1) || (fHitStrip2a[k] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }
   }

   if (fHitStrip1b[k] == 1)  {
     if (k == 1)  {
       if ((fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1) || (fHitStrip2b[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 2)  {
       if ((fHitStrip2b[k-1] == 1) ||(fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1) || (fHitStrip2b[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if ((k >= 3) && (k <= 1014))  {
       if ((fHitStrip2b[k-2] == 1) || (fHitStrip2b[k-1] == 1) || (fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1) || (fHitStrip2b[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1015)  {
       if ((fHitStrip2b[k-2] == 1) || (fHitStrip2b[k-1] == 1) || (fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1016)  {
       if ((fHitStrip2b[k-2] == 1) || (fHitStrip2b[k-1] == 1) || (fHitStrip2b[k] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }
   }
 }

 if (fHitMatchFoundInThisEvent == 1)  {
   fRunAction->AddEvtWithHitMatches();
 }

 //G4cout 
    //<< "\n totalNbHitStrips = "
    //<< totalNbHitStrips
    //<< G4endl;

 //charge weights and entrance positions of primaries in strip sensors per event
 for (G4int k = 0; k < 1017; k++) {
   fCS1a = fChargeStrip1a[k];
   //G4cout 
     //<< "\n fCS1a = "
     //<< fCS1a;
   fCS1b = fChargeStrip1b[k];
   //G4cout 
     //<< "\n fCS1b = "
     //<< fCS1b;
   fCS2a = fChargeStrip2a[k];
   //G4cout 
     //<< "\n fCS2a = "
     //<< fCS2a;
   fCS2b = fChargeStrip2b[k];
   //G4cout 
     //<< "\n fCS2b = "
     //<< fCS2b;

  if (fCS1a >= fStripThreshold)  {
      fWeightStrip1a[k] = fStripThreshold/fChargeStrip1;

      fStripCenterNRX = Strip1Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = -Dist/2 - Det1SizeZ;
      /*fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);*/

      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY;
      fStripCenterNTZ = fStripCenterNRZ;

      fStripCenterZ = fStripCenterNTZ + posZDUT;

      fPointDet1ent = fPointDet1ent + G4ThreeVector((fStripCenterX*fStripThreshold/fChargeStrip1), (fStripCenterY*fStripThreshold/fChargeStrip1), (fStripCenterZ*fStripThreshold/fChargeStrip1));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS1a = "
         //<< fCS1a;
      //G4cout 
         //<< "\n fChargeStrip1 = "
         //<< fChargeStrip1;
      //G4cout 
         //<< "\n fWeightStrip1a[k] = "
         //<< fWeightStrip1a[k];

      //G4cout 
         //<< "\n fPointDet1ent = "
         //<< fPointDet1ent;
   }
   if (fCS1b >= fStripThreshold)  {
      fWeightStrip1b[k] = fStripThreshold/fChargeStrip1;

      fStripCenterNRX = -Strip1Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = -Dist/2 - Det1SizeZ;
      /*fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);*/

      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY;
      fStripCenterNTZ = fStripCenterNRZ;

      fStripCenterZ = fStripCenterNTZ + posZDUT;

      fPointDet1ent = fPointDet1ent + G4ThreeVector((fStripCenterX*fStripThreshold/fChargeStrip1), (fStripCenterY*fStripThreshold/fChargeStrip1), (fStripCenterZ*fStripThreshold/fChargeStrip1));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS1b = "
         //<< fCS1b;
      //G4cout 
         //<< "\n fChargeStrip1 = "
         //<< fChargeStrip1;
      //G4cout 
         //<< "\n fWeightStrip1b[k] = "
         //<< fWeightStrip1b[k];

      //G4cout 
         //<< "\n fPointDet1ent = "
         //<< fPointDet1ent;
   }
   if (fCS2a >= fStripThreshold)  {
      fWeightStrip2a[k] = fStripThreshold/fChargeStrip2;

      fStripCenterNRX = Strip2Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = Dist/2;
      /*fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);*/

      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY;
      fStripCenterNTZ = fStripCenterNRZ;

      fStripCenterZ = fStripCenterNTZ + posZDUT;

      fPointDet2ent = fPointDet2ent + G4ThreeVector((fStripCenterX*fStripThreshold/fChargeStrip2), (fStripCenterY*fStripThreshold/fChargeStrip2), (fStripCenterZ*fStripThreshold/fChargeStrip2));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS2a = "
         //<< fCS2a;
      //G4cout 
         //<< "\n fChargeStrip2 = "
         //<< fChargeStrip2;
      //G4cout 
         //<< "\n fWeightStrip2a[k] = "
         //<< fWeightStrip2a[k];

      //G4cout 
         //<< "\n fPointDet2ent = "
         //<< fPointDet2ent;
   }
   if (fCS2b >= fStripThreshold)  {
      fWeightStrip2b[k] = fStripThreshold/fChargeStrip2;

      fStripCenterNRX = -Strip2Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = Dist/2;
      /*fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);*/

      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY;
      fStripCenterNTZ = fStripCenterNRZ;

      fStripCenterZ = fStripCenterNTZ + posZDUT;

      fPointDet2ent = fPointDet2ent + G4ThreeVector((fStripCenterX*fStripThreshold/fChargeStrip2), (fStripCenterY*fStripThreshold/fChargeStrip2), (fStripCenterZ*fStripThreshold/fChargeStrip2));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS2b = "
         //<< fCS2b;
      //G4cout 
         //<< "\n fChargeStrip2 = "
         //<< fChargeStrip2;
      //G4cout 
         //<< "\n fWeightStrip2b[k] = "
         //<< fWeightStrip2b[k];

      //G4cout 
         //<< "\n fPointDet2ent = "
         //<< fPointDet2ent;
   }
 }

 if ((fHitSensor1 == 1) && (fHitSensor2 == 1))  {
      fRunAction->AddNbCoinc();
 }

 // fill ntuple and charge per pixel histograms
 G4int kec = 0;
 for (G4int ke1 = 0; ke1 < 3; ke1++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
 	for (G4int ke3 = 0; ke3 < 81; ke3++) {
 	   for (G4int ke4 = 0; ke4 < 53; ke4++) {
		//analysisManager->FillNtupleDColumn(kec+8*1017, fEnergyPixel[ke1][ke2][ke3][ke4]);

   		//charge collected per pixel in number of electrons [e]
   		fChargePixel[ke1][ke2][ke3][ke4] = (fEnergyPixel[ke1][ke2][ke3][ke4])/(3.67*eV);
   		fChargePixelMod[ke1][ke2][ke3][ke4] = (fEnergyPixel[ke1][ke2][ke3][ke4])/(3.67*eV);

   		fChargePix = fChargePixel[ke1][ke2][ke3][ke4];
   		//G4cout 
     		  //<< "\n fChargePix = "
     		  //<< fChargePix;
     		  //<< " e" 
                  //<< G4endl;
 
   		if (fChargePix >= fPixThreshold)  {
      		   fRunAction->AddNbHitsPixel(ke1, ke2, ke3, ke4);
                   //fHitsMultiplicityPix[ke1] ++;
      		   //G4cout 
        	      //<< "\n NbHitsPixel advanced"
        	      //<< G4endl;
		   kec++;
   		}
 	   }
 	}
    }
 }

 //artificial charge sharing
 for (G4int ke1 = 0; ke1 < 3; ke1++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
 	for (G4int ke3 = 0; ke3 < 81; ke3++) {
 	   for (G4int ke4 = 0; ke4 < 53; ke4++) {
                fChargeQuad = (fChargePixel[ke1][ke2][ke3][ke4])/4;
		fRandPixel1 = G4RandGauss::shoot(fChargeQuad,fPixSigma);
		fRandPixel2 = G4RandGauss::shoot(fChargeQuad,fPixSigma);
		fRandPixel3 = G4RandGauss::shoot(fChargeQuad,fPixSigma);
		fRandPixel4 = G4RandGauss::shoot(fChargeQuad,fPixSigma);
		/*G4cout
		   << "\n fRandPixel1 = "
		   << fRandPixel1
		   << "\n fRandPixel2 = "
		   << fRandPixel2
		   << "\n fRandPixel3 = "
		   << fRandPixel3
		   << "\n fRandPixel4 = "
		   << fRandPixel4
		   << G4endl;*/

		iSecretPixel1 = rand() % 99;
		iSecretPixel2 = rand() % 99;
		iSecretPixel3 = rand() % 99;
		iSecretPixel4 = rand() % 99;
                if ((ke3 > 1) && (iSecretPixel1 < iRandLim))  {
		   fChargePixelMod[ke1][ke2][ke3-1][ke4] += fChargeQuad + fRandPixel1;
                }
                if ((ke3 < 80) && (iSecretPixel2 < iRandLim))  {
		   fChargePixelMod[ke1][ke2][ke3+1][ke4] += fChargeQuad + fRandPixel2;
                }
                if ((ke4 > 1) && (iSecretPixel3 < iRandLim))  {
		   fChargePixelMod[ke1][ke2][ke3][ke4-1] += fChargeQuad + fRandPixel3;
                }
                if ((ke4 < 52) && (iSecretPixel4 < iRandLim))  {
		   fChargePixelMod[ke1][ke2][ke3][ke4+1] += fChargeQuad + fRandPixel4;
                }
 	   }
 	}
    }
 }


 for (G4int ke1 = 0; ke1 < 3; ke1++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
 	for (G4int ke3 = 0; ke3 < 81; ke3++) {
 	   for (G4int ke4 = 0; ke4 < 53; ke4++) {
                fChargePixMod = fChargePixelMod[ke1][ke2][ke3][ke4];
   		if (fChargePixMod >= fPixThreshold)  {
                   fHitsMultiplicityPix[ke1] ++;
   		}
 	   }
 	}
    }
 }

 analysisManager->AddNtupleRow();

 if (fHitsMultiplicity1b > 0)  {
   analysisManager->FillH1(189, fHitsMultiplicity1b);
   analysisManager->FillH1(241, fHitsMultiplicity1b);
   //G4cout 
      //<< "\n hit multiplicity of sensor 1: "
      //<< fHitsMultiplicity1b
      //<< G4endl;
 }
 if (fHitsMultiplicity2b > 0)  {
   analysisManager->FillH1(190, fHitsMultiplicity2b);
   analysisManager->FillH1(242, fHitsMultiplicity2b);
   //G4cout 
      //<< "\n hit multiplicity of sensor 2: "
      //<< fHitsMultiplicity2b
      //<< G4endl;
 }

 if (fHitsMultiplicityPix[1] > 0)  {
   analysisManager->FillH1(197, fHitsMultiplicityPix[1]);
   analysisManager->FillH1(239, fHitsMultiplicityPix[1]);
   //G4cout 
      //<< "\n hit multiplicity of BPIX module 1: "
      //<< fHitsMultiplicityPix[1]
      //<< G4endl;
 }
 if (fHitsMultiplicityPix[2] > 0)  {
   analysisManager->FillH1(198, fHitsMultiplicityPix[2]);
   analysisManager->FillH1(240, fHitsMultiplicityPix[2]);
   //G4cout 
      //<< "\n hit multiplicity of BPIX module 2: "
      //<< fHitsMultiplicityPix[2]
      //<< G4endl;
 }

 fCharge1 = fEnergyDeposit1/(3.67*eV);
 fCharge2 = fEnergyDeposit2/(3.67*eV);
 analysisManager->FillH1(155, fCharge1);
 analysisManager->FillH1(156, fCharge2);

 for (G4int ke2 = 0; ke2 < 17; ke2++) {
   for (G4int ke3 = 0; ke3 < 81; ke3++) {
     for (G4int ke4 = 0; ke4 < 53; ke4++) {
       //charge collected per pixel in number of electrons [e]
       fChargePixel[1][ke2][ke3][ke4] = (fEnergyPixel[1][ke2][ke3][ke4])/(3.67*eV);
       fChargePixel[2][ke2][ke3][ke4] = (fEnergyPixel[2][ke2][ke3][ke4])/(3.67*eV);

       fChargePix1 = fChargePixel[1][ke2][ke3][ke4];
       fChargePix2 = fChargePixel[2][ke2][ke3][ke4];
 
       if (fChargePix1 >= fPixThreshold)  {
	 fChargePixel1 = fChargePixel1 + fChargePix1;
 	 //G4cout 
    	    //<< "\n Pixel hit: Module: 1 ROC: "
    	    //<< ke2
    	    //<< " row: "
    	    //<< ke3
    	    //<< " col: "
    	    //<< ke4
    	    //<< G4endl;
       }
       if (fChargePix2 >= fPixThreshold)  {
	 fChargePixel2 = fChargePixel2 + fChargePix2;
 	 //G4cout 
    	    //<< "\n Pixel hit: Module: 2 ROC: "
    	    //<< ke2
    	    //<< " row: "
    	    //<< ke3
    	    //<< " col: "
    	    //<< ke4
    	    //<< G4endl;
       }
     }
   }
 }


 //charge weights and entrance positions of primaries in BPIX modules per event
 for (G4int ke2 = 0; ke2 < 17; ke2++) {
   for (G4int ke3 = 0; ke3 < 81; ke3++) {
     for (G4int ke4 = 0; ke4 < 53; ke4++) {
       fChargePix1 = fChargePixel[1][ke2][ke3][ke4];
       fChargePix2 = fChargePixel[2][ke2][ke3][ke4];
 
       if (fChargePix1 >= fPixThreshold)  {
	 fHitPixelDet1 = 1;
         fWeightPixel1[ke2][ke3][ke4] = fChargePix1/fChargePixel1;

	 if ((ke2 == 1) || (ke2 == 9))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 3*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 3*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 3*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 2) || (ke2 == 10))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 2*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 2*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 2*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 3) || (ke2 == 11))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 1*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 1*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 1*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 4) || (ke2 == 12))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 0*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 0*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 0*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 5) || (ke2 == 13))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -0*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -0*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -0*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 6) || (ke2 == 14))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -1*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -1*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -1*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 7) || (ke2 == 15))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -2*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -2*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -2*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 8) || (ke2 == 16))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -3*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -3*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -3*54*pixelX - 53*pixelX;
	   }
	 }

	 if ((ke2>=1) && (ke2<=8))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = (81 - ke3)*pixelY + pixelY/2;
           }
	 }
	 if ((ke2>=9) && (ke2<=16))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = -80*pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = -(ke3 - 1)*pixelY - pixelY/2;
           }
	 }

	 fPixelCenterNRZ = posZBPIX12 - BPIXSizeZ/2;	 

  	 fPixelCenterNRX1 = fPixelCenterNRX;
  	 fPixelCenterNRY1 = fPixelCenterNRY;
   	 fPixelCenterNRZ1 = fPixelCenterNRZ;

  	 fPixelCenterX = fPixelCenterNRX;
  	 fPixelCenterY = fPixelCenterNRY;
  	 fPixelCenterZ = fPixelCenterNRZ;

         fPointPix1ent = fPointPix1ent + G4ThreeVector((fPixelCenterX*fChargePix1/fChargePixel1), (fPixelCenterY*fChargePix1/fChargePixel1), (fPixelCenterZ*fChargePix1/fChargePixel1));
       }

       if (fChargePix2 >= fPixThreshold)  {
	 fHitPixelDet2 = 1;
         fWeightPixel2[ke2][ke3][ke4] = fChargePix2/fChargePixel2;

	 if ((ke2 == 1) || (ke2 == 9))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 3*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 3*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 3*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 2) || (ke2 == 10))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 2*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 2*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 2*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 3) || (ke2 == 11))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 1*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 1*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 1*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 4) || (ke2 == 12))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 0*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 0*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 0*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 5) || (ke2 == 13))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -0*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -0*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -0*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 6) || (ke2 == 14))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -1*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -1*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -1*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 7) || (ke2 == 15))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -2*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -2*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -2*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 8) || (ke2 == 16))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -3*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -3*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -3*54*pixelX - 53*pixelX;
	   }
	 }

	 if ((ke2>=1) && (ke2<=8))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = (81 - ke3)*pixelY + pixelY/2;
           }
	 }
	 if ((ke2>=9) && (ke2<=16))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = -80*pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = -(ke3 - 1)*pixelY - pixelY/2;
           }
	 }

	 fPixelCenterNRZ = posZBPIX34 - BPIXSizeZ/2;	 

  	 fPixelCenterNRX1 = fPixelCenterNRX;
  	 fPixelCenterNRY1 = fPixelCenterNRY;
   	 fPixelCenterNRZ1 = fPixelCenterNRZ;

  	 fPixelCenterX = fPixelCenterNRX;
  	 fPixelCenterY = fPixelCenterNRY;
  	 fPixelCenterZ = fPixelCenterNRZ;

         fPointPix2ent = fPointPix2ent + G4ThreeVector((fPixelCenterX*fChargePix2/fChargePixel2), (fPixelCenterY*fChargePix2/fChargePixel2), (fPixelCenterZ*fChargePix2/fChargePixel2));
       }
     }
   }
 }

 for (G4int i2 = 157; i2 < 165; i2 = i2 + 1) {
   analysisManager->FillH1(i2, fChargeStrip1a[i2+348]);
 }
 for (G4int p2 = 165; p2 < 173; p2 = p2 + 1) {
   analysisManager->FillH1(p2, fChargeStrip1b[p2+340]);
 }
 for (G4int l2 = 173; l2 < 181; l2 = l2 + 1) {
   analysisManager->FillH1(l2, fChargeStrip2a[l2+332]);
 }
 for (G4int n2 = 181; n2 < 189; n2 = n2 + 1) {
   analysisManager->FillH1(n2, fChargeStrip2b[n2+324]);
 }

 ftheta = acos(fMomDir1x*fMomDir2x + fMomDir1y*fMomDir2y + fMomDir1z*fMomDir2z);
 analysisManager->FillH1(116, ftheta);


 /*//cluster occupancy
 for (G4int kea = 0; kea < 3; kea++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
      for (G4int ke3 = 0; ke3 < 81; ke3++) {// y-direction, two modules
        for (G4int ke4 = 0; ke4 < 53; ke4++) {// x-direction, one module
          fChargePix = fChargePixel[kea][ke2][ke3][ke4];
 
          if (fChargePix >= fPixThreshold)  {
	     if ((ke2 >= 9) && (ke2 <= 16))   {
		mod = 2*kea;
	        row = ke3;
	        col = (16 - ke2)*52 + 53 - ke4;
	     }
	     if ((ke2 >= 1) && (ke2 <= 8))   {
		mod = 2*kea - 1;
	        row = 80 + ke3;
	        col = (8 - ke2)*52 + 53 - ke4;
	     }
	     fClusterOccupancy[mod][row][col]++;
   	     if (kea > 0)  {
      	        analysisManager->FillH2((mod+10), col, row);
   	     }
          }
	}
     }
   }
 }*/

 //cluster occupancy
 for (G4int kea = 1; kea < 3; kea++) {
    for (G4int ke2 = 1; ke2 < 17; ke2++) {
      for (G4int ke3 = 1; ke3 < 81; ke3++) {// y-direction, two modules
        for (G4int ke4 = 1; ke4 < 53; ke4++) {// x-direction, one module
          fChargePix = fChargePixel[kea][ke2][ke3][ke4];
 
          if (fChargePix >= fPixThreshold)  {
	     /*if ((ke2 >= 9) && (ke2 <= 16))   {
		mod = 2*kea - 1;
	        row = 161 - ke3;
	        col = (16 - ke2)*52 + 53 - ke4;
	     }
	     if ((ke2 >= 1) && (ke2 <= 8))   {
		mod = 2*kea;
	        row = 81 - ke3;
	        col = (8 - ke2)*52 + 53 - ke4;
	     }*/

	     if ((ke2 >= 9) && (ke2 <= 16))   {
		mod = 2*kea - 1;
	        row = 81 - ke3;
	        col = (16 - ke2)*52 + 53 - ke4;
	     }
	     if ((ke2 >= 1) && (ke2 <= 8))   {
		mod = 2*kea - 1;
	        row = 161 - ke3;
	        col = (8 - ke2)*52 + 53 - ke4;
	     }

	     fClusterOccupancy[mod][row][col]++;
   	     if ((kea > 0) && (mod > 0))  {
      	        analysisManager->FillH2((mod+10), col, row);
      	        analysisManager->FillH1((mod+242), col);
      	        analysisManager->FillH1((mod+243), row);
   	     }
          }
	}
     }
   }
 }

 //G4cout
    //<< "\n The deflection angle is: "
    //<< ftheta
    //<< " rad"
    //<< G4endl;

 xBA = fPointDet1ent - fPointPix1ent;
 xBD = fPointDet1ent - fPointPix2ent;
 xCA = fPointDet2ent - fPointPix1ent;
 xCD = fPointDet2ent - fPointPix2ent;
 xAD = fPointPix2ent - fPointPix1ent;

 xBAD = G4ThreeVector((xBA.y())*(xBD.z()) - (xBA.z())*(xBD.y()), (xBA.z())*(xBD.x()) - (xBA.x())*(xBD.z()), (xBA.x())*(xBD.y()) - (xBA.y())*(xBD.x()));
 xCAD = G4ThreeVector((xCA.y())*(xCD.z()) - (xCA.z())*(xCD.y()), (xCA.z())*(xCD.x()) - (xCA.x())*(xCD.z()), (xCA.x())*(xCD.y()) - (xCA.y())*(xCD.x()));

 dBAD = sqrt(sqr(xBAD.x()) + sqr(xBAD.y()) + sqr(xBAD.z()));
 dCAD = sqrt(sqr(xCAD.x()) + sqr(xCAD.y()) + sqr(xCAD.z()));
 dAD = sqrt(sqr(xAD.x()) + sqr(xAD.y()) + sqr(xAD.z()));

 dB = dBAD/dAD;
 dC = dCAD/dAD;

 //G4cout 
    //<< "\n Distance between B and AD: " 
    //<< G4BestUnit(dB, "Length")
    //<< G4endl;
 if (dB == dB)   {
    analysisManager->FillH1(117, dB);
 }

 //G4cout 
    //<< "\n Distance between C and AD: " 
    //<< G4BestUnit(dC, "Length")
    //<< G4endl;
 if (dC == dC)   {
    analysisManager->FillH1(118, dC);
 }

 B1 = G4ThreeVector((xAD.x())*(fPointDet1ent.z() - fPointPix2ent.z())/(xAD.z()) + fPointPix2ent.x(), (xAD.y())*(fPointDet1ent.z() - fPointPix2ent.z())/(xAD.z()) + fPointPix2ent.y(), fPointDet1ent.z());
 B1Bx = fPointDet1ent.x() - B1.x();
 B1By = fPointDet1ent.y() - B1.y();


 if (B1Bx == B1Bx)  {
   analysisManager->FillH1(119, B1Bx);

 //G4cout 
    //<< "\n fPointDet1ent.x(): " 
    //<< G4BestUnit(fPointDet1ent.x(), "Length")
    //<< "\n B1.x(): " 
    //<< G4BestUnit(B1.x(), "Length")
    //<< G4endl;

 }
 if (B1By == B1By)  {
   analysisManager->FillH1(120, B1By);
   if ((fPointDet1ent.y() > -2*StripPitch/3) && (fPointDet1ent.y() < -StripPitch/90))   {
      analysisManager->FillH1(209, B1By);
   }
   if ((fPointDet1ent.y() > -StripPitch/90) && (fPointDet1ent.y() < StripPitch/90))   {
      analysisManager->FillH1(210, B1By);
   }
   if ((fPointDet1ent.y() > StripPitch/90) && (fPointDet1ent.y() < 2*StripPitch/3))   {
      analysisManager->FillH1(211, B1By);
   }
 }

 C1 = G4ThreeVector((xAD.x())*(fPointDet2ent.z() - fPointPix1ent.z())/(xAD.z()) + fPointPix1ent.x(), (xAD.y())*(fPointDet2ent.z() - fPointPix1ent.z())/(xAD.z()) + fPointPix1ent.y(), fPointDet2ent.z());
 C1Cx = fPointDet2ent.x() - C1.x();
 C1Cy = fPointDet2ent.y() - C1.y();
 //G4cout 
    //<< "\n C'Cx: " 
    //<< G4BestUnit(C1Cx, "Length")
    //<< "\n C'Cy: " 
    //<< G4BestUnit(C1Cy, "Length")
    //<< G4endl;

 if (C1Cx == C1Cx)  {
   analysisManager->FillH1(121, C1Cx);
 }
 if (C1Cy == C1Cy)  {
   analysisManager->FillH1(122, C1Cy);
 }

 /*G4cout
    << "\n fPointDet1ent = "
    << G4BestUnit(fPointDet1ent, "Length")
    << "\n fPointDet2ent = "
    << G4BestUnit(fPointDet2ent, "Length")
    << "\n fPointPix1ent = "
    << G4BestUnit(fPointPix1ent, "Length")
    << "\n fPointPix2ent = "
    << G4BestUnit(fPointPix2ent, "Length")
    << "\n B1 = "
    << G4BestUnit(B1, "Length")
    << "\n C1 = "
    << G4BestUnit(C1, "Length")
    << "\n B1Bx = "
    << G4BestUnit(B1Bx, "Length")
    << "\n B1By = "
    << G4BestUnit(B1By, "Length")
    << "\n C1Cx = "
    << G4BestUnit(C1Cx, "Length")
    << "\n C1Cy = "
    << G4BestUnit(C1Cy, "Length")
    << G4endl;*/

 // check for stubs in 2S sensors
 if (B1By <= stubLim)  {
   for (G4int k = 0; k < 1017; k++) {
     if (fHitStrip1a[k] == 1)  {
       if (k == 1)  {
         if ((fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1) || (fHitStrip2a[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 2)  {
       if ((fHitStrip2a[k-1] == 1) ||(fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1) || (fHitStrip2a[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if ((k >= 3) && (k <= 1014))  {
       if ((fHitStrip2a[k-2] == 1) || (fHitStrip2a[k-1] == 1) || (fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1) || (fHitStrip2a[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1015)  {
       if ((fHitStrip2a[k-2] == 1) || (fHitStrip2a[k-1] == 1) || (fHitStrip2a[k] == 1) || (fHitStrip2a[k+1] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1016)  {
       if ((fHitStrip2a[k-2] == 1) || (fHitStrip2a[k-1] == 1) || (fHitStrip2a[k] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }
   }

   if (fHitStrip1b[k] == 1)  {
     if (k == 1)  {
       if ((fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1) || (fHitStrip2b[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 2)  {
       if ((fHitStrip2b[k-1] == 1) ||(fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1) || (fHitStrip2b[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if ((k >= 3) && (k <= 1014))  {
       if ((fHitStrip2b[k-2] == 1) || (fHitStrip2b[k-1] == 1) || (fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1) || (fHitStrip2b[k+2] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1015)  {
       if ((fHitStrip2b[k-2] == 1) || (fHitStrip2b[k-1] == 1) || (fHitStrip2b[k] == 1) || (fHitStrip2b[k+1] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }

     if (k == 1016)  {
       if ((fHitStrip2b[k-2] == 1) || (fHitStrip2b[k-1] == 1) || (fHitStrip2b[k] == 1)) {
          fHitMatchFoundInThisEvent = 1;
          fRunAction->AddHitMatches();
	  /*G4cout
	     << "\n A hit match has been found."
	     << G4endl;*/
       }
     }
   }
 }
 }


 //G4cout
    //<< "\n A.x = "
    //<< G4BestUnit(fPointPix1ent.x(), "Length")
    //<< "\n A.y = "
    //<< G4BestUnit(fPointPix1ent.y(), "Length")
    //<< "\n A.z = "
    //<< G4BestUnit(fPointPix1ent.z(), "Length")
    //<< "\n B.x = "
    //<< G4BestUnit(fPointDet1ent.x(), "Length")
    //<< "\n B.y = "
    //<< G4BestUnit(fPointDet1ent.y(), "Length")
    //<< "\n B.z = "
    //<< G4BestUnit(fPointDet1ent.z(), "Length")
    //<< "\n B'.x = "
    //<< G4BestUnit(B1.x(), "Length")
    //<< "\n B'.y = "
    //<< G4BestUnit(B1.y(), "Length")
    //<< "\n B'.z = "
    //<< G4BestUnit(B1.z(), "Length")
    //<< "\n C.x = "
    //<< G4BestUnit(fPointDet2ent.x(), "Length")
    //<< "\n C.y = "
    //<< G4BestUnit(fPointDet2ent.y(), "Length")
    //<< "\n C.z = "
    //<< G4BestUnit(fPointDet2ent.z(), "Length")
    //<< "\n C'.x = "
    //<< G4BestUnit(C1.x(), "Length")
    //<< "\n C'.y = "
    //<< G4BestUnit(C1.y(), "Length")
    //<< "\n C'.z = "
    //<< G4BestUnit(C1.z(), "Length")
    //<< "\n D.x = "
    //<< G4BestUnit(fPointPix2ent.x(), "Length")
    //<< "\n D.y = "
    //<< G4BestUnit(fPointPix2ent.y(), "Length")
    //<< "\n D.z = "
    //<< G4BestUnit(fPointPix2ent.z(), "Length")
    //<< "\n \n fPointPix1entReal.x = "
    //<< G4BestUnit(fPointPix1entReal.x(), "Length")
    //<< "\n fPointPix1entReal.y = "
    //<< G4BestUnit(fPointPix1entReal.y(), "Length")
    //<< "\n fPointPix1entReal.z = "
    //<< G4BestUnit(fPointPix1entReal.z(), "Length")
    //<< "\n fPointDet1entReal.x = "
    //<< G4BestUnit(fPointDet1entReal.x(), "Length")
    //<< "\n fPointDet1entReal.y = "
    //<< G4BestUnit(fPointDet1entReal.y(), "Length")
    //<< "\n fPointDet1entReal.z = "
    //<< G4BestUnit(fPointDet1entReal.z(), "Length")
    //<< "\n fPointDet2entReal.x = "
    //<< G4BestUnit(fPointDet2entReal.x(), "Length")
    //<< "\n fPointDet2entReal.y = "
    //<< G4BestUnit(fPointDet2entReal.y(), "Length")
    //<< "\n fPointDet2entReal.z = "
    //<< G4BestUnit(fPointDet2entReal.z(), "Length")
    //<< "\n fPointPix2entReal.x = "
    //<< G4BestUnit(fPointPix2entReal.x(), "Length")
    //<< "\n fPointPix2entReal.y = "
    //<< G4BestUnit(fPointPix2entReal.y(), "Length")
    //<< "\n fPointPix2entReal.z = "
    //<< G4BestUnit(fPointPix2entReal.z(), "Length")
    //<< G4endl;


 fDiff1x = fPointDet1ent.x() - fPointDet1entReal.x();
 fDiff1y = fPointDet1ent.y() - fPointDet1entReal.y();
 fDiff2x = fPointDet2ent.x() - fPointDet2entReal.x();
 fDiff2y = fPointDet2ent.y() - fPointDet2entReal.y();

 if (fHitSensor1 == 1)   {
    analysisManager->FillH1(191, fPointDet1ent.x());
    analysisManager->FillH1(192, fPointDet1ent.y());
    analysisManager->FillH1(193, fPointDet1ent.z());
    analysisManager->FillH1(199, fDiff1x);
    analysisManager->FillH1(200, fDiff1y);

    //analysisManager->FillH2(1, fPointDet1ent.y(), fPointDet1entReal.y());
 }
 if (fHitSensor2 == 1)   {
    analysisManager->FillH1(194, fPointDet2ent.x());
    analysisManager->FillH1(195, fPointDet2ent.y());
    analysisManager->FillH1(196, fPointDet2ent.z());
    analysisManager->FillH1(201, fDiff2x);
    analysisManager->FillH1(202, fDiff2y);

    //analysisManager->FillH2(2, fPointDet2ent.y(), fPointDet2entReal.y());
 }
 if (fHitPixelDet1 == 1)   {
    analysisManager->FillH1(203, fPointPix1ent.x());
    analysisManager->FillH1(204, fPointPix1ent.y());
    analysisManager->FillH1(205, fPointPix1ent.z());
 }
 if (fHitPixelDet2 == 1)   {
    analysisManager->FillH1(206, fPointPix2ent.x());
    analysisManager->FillH1(207, fPointPix2ent.y());
    analysisManager->FillH1(208, fPointPix2ent.z());
 }

 analysisManager->FillH2(1, fPoint20Scint.x(), fPoint20Scint.y()-8.0*mm);
 analysisManager->FillH2(2, fPointExitAl.x(), fPointExitAl.y()-8.0*mm);
 analysisManager->FillH2(3, fPointExitScint1.x(), fPointExitScint1.y()-8.0*mm);
 analysisManager->FillH2(4, fPointExitPixel12.x(), fPointExitPixel12.y()-8.0*mm);
 analysisManager->FillH2(5, fPointExitSensor1.x(), fPointExitSensor1.y()-8.0*mm);
 analysisManager->FillH2(6, fPointExitSensor2.x(), fPointExitSensor2.y()-8.0*mm);
 analysisManager->FillH2(7, fPointExitPixel34.x(), fPointExitPixel34.y()-8.0*mm);
 analysisManager->FillH2(8, fPointExitScint2.x(), fPointExitScint2.y()-8.0*mm);
 analysisManager->FillH2(9, fEnergyDepositSc2, fEnergyEntScintillator2);

 //G4cout 
    //<< "\n fPointDet1ent.x() = " 
    //<< G4BestUnit(fPointDet1ent.x(), "Length")
    //<< "\n fPointDet1entReal.x() = " 
    //<< G4BestUnit(fPointDet1entReal.x(), "Length")
    //<< G4endl;

 //G4cout 
     //<< "\n Last track length of secondary created in detector calculated = " 
     //<< G4BestUnit(fTrack2, "Length")
     //<< "\n End of Event. Secondary track length from detector secondaries = " 
     //<< G4BestUnit(fSecondaryDetTrackLength, "Length")
     //<< G4endl;

 //Visualize event if there is a track longer than 1 cm
 //if (fTrack2 > 1.0*cm)  {
     //G4cout 
         //<< "\n fTrack2 = " 
         //<< G4BestUnit(fTrack2, "Length")
         //<< G4endl;
     //G4EventManager* evMan = G4EventManager::GetEventManager();
     //evMan->KeepTheCurrentEvent();
 //}

 //G4cout 
    //<< "\n End of event! "
    //<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

