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
/// \file electromagnetic/TestEm18/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*, RunAction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void   EndOfEventAction(const G4Event*);

    
    void AddEnergyDeposit1(G4double edep)                     {fEnergyDeposit1  += edep;};
    void AddEnergyDeposit2(G4double edep)                     {fEnergyDeposit2  += edep;};
    void AddEnergyDepositL1(G4double edep)                    {fEnergyDepositL1  += edep;};
    void AddEnergyDepositL2(G4double edep)                    {fEnergyDepositL2  += edep;};
    void AddEnergyDepositROCL1(G4double edep)                 {fEnergyDepositROCL1  += edep;};
    void AddEnergyDepositROCL2(G4double edep)                 {fEnergyDepositROCL2  += edep;};
    void AddEnergyDepositSc1(G4double edep)                   {fEnergyDepositSc1  += edep;};
    void AddEnergyDepositSc2(G4double edep)                   {fEnergyDepositSc2  += edep;};
    void AddEnergyDepositAir(G4double edep)                   {fEnergyDepositAir  += edep;};
    void AddEnergyDepositAl(G4double edep)                    {fEnergyDepositAl  += edep;};

    void AddEnergyEntSensor1(G4double edep)                   {fEnergyEntSensor1  += edep;};
    void AddEnergyExSensor1(G4double edep)                    {fEnergyExSensor1  += edep;};
    void AddEnergyEntSensor2(G4double edep)                   {fEnergyEntSensor2  += edep;};
    void AddEnergyExSensor2(G4double edep)                    {fEnergyExSensor2  += edep;};
    void AddEnergyEntPixelLayer1(G4double edep)               {fEnergyEntPixelLayer1  += edep;};
    void AddEnergyExPixelLayer1(G4double edep)                {fEnergyExPixelLayer1  += edep;};
    void AddEnergyEntPixelLayer2(G4double edep)               {fEnergyEntPixelLayer2  += edep;};
    void AddEnergyExPixelLayer2(G4double edep)                {fEnergyExPixelLayer2  += edep;};
    void AddEnergyEntScintillator1(G4double edep)             {fEnergyEntScintillator1  += edep;};
    void AddEnergyExScintillator1(G4double edep)              {fEnergyExScintillator1  += edep;};
    void AddEnergyEntScintillator2(G4double edep)             {fEnergyEntScintillator2  += edep;};
    void AddEnergyExScintillator2(G4double edep)              {fEnergyExScintillator2  += edep;};
    void AddEnergyExAl(G4double edep)                         {fEnergyExAl  += edep;};

    void AddPrimaryTrackLength(G4double track)                {fPrimaryTrackLength  += track;};
    void AddSecondary1(G4double ekin)                         {fEnergySecondary1  += ekin;};
    void AddSecondary2(G4double ekin)                         {fEnergySecondary2  += ekin;};
    void AddSecondaryxPolarization(G4double spolarization)    {fSecondaryxPolarization  += spolarization;};
    void AddSecondaryyPolarization(G4double spolarization)    {fSecondaryyPolarization  += spolarization;};
    void AddSecondaryzPolarization(G4double spolarization)    {fSecondaryzPolarization  += spolarization;};
    void AddSecondaryTrackLength(G4double track)              {fSecondaryTrackLength  += track;};
    void AddSecondaryDetTrackLength(G4double track)           {fSecondaryDetTrackLength  += track;};
    void AddTertiary1(G4double ekin)                          {fEnergyTertiary1  += ekin;};
    void AddTertiary2(G4double ekin)                          {fEnergyTertiary2  += ekin;};
    void AddTertiaryxPolarization(G4double spolarization)     {fTertiaryxPolarization  += spolarization;};
    void AddTertiaryyPolarization(G4double spolarization)     {fTertiaryyPolarization  += spolarization;};
    void AddTertiaryzPolarization(G4double spolarization)     {fTertiaryzPolarization  += spolarization;};
    void AddTertiaryTrackLength(G4double track)               {fTertiaryTrackLength  += track;};
    void TrackCheck(G4double trackcheck)                      {fTrack1 = fTrack2; fTrack2 = trackcheck; 
                                                               if (fTrack1 > fTrack2) {
                                                                   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
                                                                   analysisManager->FillH1(31, fTrack1);
                                                                   //G4cout 
     									//<< "\n Last track length of secondary created in detector calculated = " 
     									//<< G4BestUnit(fTrack1, "Length")
     									//<< G4endl;
                                                               }
							      };
    void AddEnergyStrip1a(G4double edep, G4int i)  {fEnergyStrip1a[i] += edep;};
    void AddEnergyStrip1b(G4double edep, G4int i)  {fEnergyStrip1b[i] += edep;};
    void AddEnergyStrip2a(G4double edep, G4int i)  {fEnergyStrip2a[i] += edep;};
    void AddEnergyStrip2b(G4double edep, G4int i)  {fEnergyStrip2b[i] += edep;};

    void AddEnergyPixel(G4double edep, G4int i, G4int j, G4int k, G4int l)  {fEnergyPixel[i][j][k][l] += edep;};

    void AddMomentumDirection1(G4double dirx, G4double diry, G4double dirz)  {fMomDir1x += dirx; fMomDir1y += diry; fMomDir1z += dirz;};
    void AddMomentumDirection2(G4double dirx, G4double diry, G4double dirz)  {fMomDir2x += dirx; fMomDir2y += diry; fMomDir2z += dirz;};

    void AddPointDet1entReal(G4ThreeVector vec)  {fPointDet1entReal += vec;};
    void AddPointDet2entReal(G4ThreeVector vec)  {fPointDet2entReal += vec;};

    void AddPointPix1entReal(G4ThreeVector vec)  {fPointPix1entReal += vec;};
    void AddPointPix2entReal(G4ThreeVector vec)  {fPointPix2entReal += vec;};

    void AddPointAfterScint2(G4ThreeVector vec)  {fPoint20Scint += vec;};

    void AddPointExitAl(G4ThreeVector vec)  {fPointExitAl += vec;};
    void AddPointExitScint1(G4ThreeVector vec)  {fPointExitScint1 += vec;};
    void AddPointExitScint2(G4ThreeVector vec)  {fPointExitScint2 += vec;};
    void AddPointExitPixel12(G4ThreeVector vec)  {fPointExitPixel12 += vec;};
    void AddPointExitPixel34(G4ThreeVector vec)  {fPointExitPixel34 += vec;};
    void AddPointExitSensor1(G4ThreeVector vec)  {fPointExitSensor1 += vec;};
    void AddPointExitSensor2(G4ThreeVector vec)  {fPointExitSensor2 += vec;};

    void AddTimeEntranceScint1(G4double time)  {fTimeEntranceScint1 += time;};
    void AddTimeEntranceScint2(G4double time)  {fTimeEntranceScint2 += time;};

    void AddNb2SElectrons()  {fNb2SElectrons++;
			      /*G4cout
				<< "\n fNb2SElectrons++"
				<< fNb2SElectrons
				<< G4endl;*/};


        
  private:
    DetectorConstruction* fDetectorconstruction;
    RunAction*    fRunAction;
    
    G4double      fEnergyDeposit1;
    G4double      fEnergyDeposit2;
    G4double      fEnergyDepositL1;
    G4double      fEnergyDepositL2;
    G4double      fEnergyDepositROCL1;
    G4double      fEnergyDepositROCL2;
    G4double      fEnergyDepositSc1;
    G4double      fEnergyDepositSc2;
    G4double      fEnergyDepositAir;
    G4double      fEnergyDepositAl;

    G4double      fEnergyEntSensor1;
    G4double      fEnergyExSensor1;
    G4double      fEnergyEntSensor2;
    G4double      fEnergyExSensor2;
    G4double      fEnergyEntPixelLayer1;
    G4double      fEnergyExPixelLayer1;
    G4double      fEnergyEntPixelLayer2;
    G4double      fEnergyExPixelLayer2;
    G4double      fEnergyEntScintillator1;
    G4double      fEnergyExScintillator1;
    G4double      fEnergyEntScintillator2;
    G4double      fEnergyExScintillator2;
    G4double      fEnergyExAl;

    G4double      fPrimaryTrackLength;
    G4double      fEnergySecondary1;
    G4double      fEnergySecondary2;       
    G4double      fSecondaryxPolarization;  
    G4double      fSecondaryyPolarization;  
    G4double      fSecondaryzPolarization; 
    G4double      fSecondaryTrackLength; 
    G4double      fSecondaryDetTrackLength;
    G4double      fEnergyTertiary1;
    G4double      fEnergyTertiary2;   
    G4double      fTertiaryxPolarization;   
    G4double      fTertiaryyPolarization; 
    G4double      fTertiaryzPolarization; 
    G4double      fTertiaryTrackLength;
    G4double      fEnergyStrip1a[1017];
    G4double      fEnergyStrip1b[1017];
    G4double      fEnergyStrip2a[1017];
    G4double      fEnergyStrip2b[1017];
    G4double      fWeightStrip1a[1017];
    G4double      fWeightStrip1b[1017];
    G4double      fWeightStrip2a[1017];
    G4double      fWeightStrip2b[1017];
    G4int         fHitStrip1a[1017];
    G4int         fHitStrip1b[1017];
    G4int         fHitStrip2a[1017];
    G4int         fHitStrip2b[1017];
    G4double      fMomDir1x, fMomDir1y, fMomDir1z;
    G4double      fMomDir2x, fMomDir2y, fMomDir2z;
    G4double      fEnergyPixel[3][17][81][53];
    G4double      fWeightPixel1[17][81][53];
    G4double      fWeightPixel2[17][81][53];

    G4double 	  fStripCenterNRX, fStripCenterNRY, fStripCenterNRZ, fStripCenterNTZ; //Strip center position if the DUT wasn't rotated or translated
    G4double 	  fStripCenterX, fStripCenterY, fStripCenterZ;

    G4double 	  fPixelCenterNRX, fPixelCenterNRY, fPixelCenterNRZ, fPixelCenterNTZ, fPixelCenterNRX1, fPixelCenterNRY1, fPixelCenterNRZ1; //Pixel center position if the DUT wasn't rotated
    G4double 	  fPixelCenterX, fPixelCenterY, fPixelCenterZ;

    G4double	  fDiff1x, fDiff1y, fDiff2x, fDiff2y;

    G4double      Det1SizeZ, Det2SizeZ, Dist, Strip1Depth, Strip1Length, Strip2Depth, Strip2Length, StripDist, StripWidth, StripPitch, posEndArm1, posBeginningArm2, BPIXSizeZ, pixelDepth, pixelX, pixelY, ROChor, ROCvert, XangleDUT, XangleBPIX, YangleBPIX, ElField1, ElField2;

    G4double      posZAl, posZScint1, posZScint2, posZBPIX12, posZBPIX34, posZ1, posZ2, posZDUT;

    G4int	  fRandStrip1, fRandStrip2, fRandStrip3, fRandStrip4, fRandStrip5, fRandStrip6, fRandStrip7, fRandStrip8, fRandPixel1, fRandPixel2, fRandPixel3, fRandPixel4;

    G4int 	  totalNbStrips, stripNo, totalNbHitStrips;
    G4int         fCharge1, fCharge2, fCS1a, fCS1b, fCS2a, fCS2b, fChargePix, fChargePixMod, fChargeStrip1, fChargeStrip2, fChargePix1, fChargePix2, fChargePixel1, fChargePixel2, fChargeQuad;
    G4int         fHitSensor1, fHitSensor2;
    G4int	  fHitPixelDet1, fHitPixelDet2;

    G4int         fNbHitsStrip1a[1017];
    G4int         fNbHitsStrip1b[1017];
    G4int         fNbHitsStrip2a[1017];
    G4int         fNbHitsStrip2b[1017];
    G4int         fChargeStrip1a[1017];
    G4int         fChargeStripMod1a[1017];
    G4int         fChargeStrip1b[1017];
    G4int         fChargeStripMod1b[1017];
    G4int         fChargeStrip2a[1017];
    G4int         fChargeStripMod2a[1017];
    G4int         fChargeStrip2b[1017];
    G4int         fChargeStripMod2b[1017];
    G4int         fChargePixel[3][17][81][53];
    G4int         fChargePixelMod[3][17][81][53];
    G4int         fNbHitsPixel[3][17][81][53];
    G4int         fNb2SElectrons;

    G4int         fHitsMultiplicity1b;
    G4int         fHitsMultiplicity2b;
    G4int         fHitsMultiplicityPix[3];

    G4int         mod, row, col;
    G4int         fClusterOccupancy[5][417][161];

    G4int         fHitMatchFoundInThisEvent, fStubFoundInThisEvent;
    G4int         fPixThreshold, fStripThreshold, fPixSigma, fStripSigma;

    G4int	  iSecretStrip1, iSecretStrip2, iSecretStrip3, iSecretStrip4, iSecretStrip5, iSecretStrip6, iSecretStrip7, iSecretStrip8, iSecretPixel1, iSecretPixel2, iSecretPixel3, iSecretPixel4, iRandLim;

    G4ThreeVector fPointPix1ent, fPointPix1entReal, fPointPix1mid, fPointPix1ex, fPointDet1ent, fPointDet1entReal, fPointDet1mid, fPointDet1ex, fPointDet2ent, fPointDet2entReal, fPointDet2mid, fPointDet2ex, fPointPix2ent, fPointPix2entReal,  fPointPix2mid, fPointPix2ex, fPoint20Scint, fPointExitAl, fPointExitScint1, fPointExitScint2, fPointExitPixel12, fPointExitPixel34, fPointExitSensor1, fPointExitSensor2;

    G4double      fTimeEntranceScint1, fTimeEntranceScint2;

    //Let A = fPointPix1mid, B = fPointDet1mid, C = fPointDet2mid, D = fPointPix2mid. Let dB = distance between B and AD and dC = distance between C and AD.
    //Let B1 be the AD point with same z as B. Let C1 be the AD point with the same z as C.
    G4double dB, dC;
    G4double dBAD, dCAD, dAD;
    G4ThreeVector xBA, xBD, xCA, xCD, xAD;
    G4ThreeVector xBAD, xCAD;
    G4ThreeVector B1, C1;
    G4double B1Bx, B1By, C1Cx, C1Cy;


    //primary track's deflection angle in radians
    G4double ftheta;

    G4double      fTrack1;
    G4double      fTrack2;

    G4double      stubLim;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
