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
// $Id: SteppingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "StackingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4VTrajectory.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Gamma.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* DA, RunAction* RA, EventAction* EA, StackingAction* SA)
:G4UserSteppingAction(), fDetectorconstruction(DA), fRunaction(RA), fEventaction(EA), fStackingaction(SA)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 //get volume of the current step
  G4VPhysicalVolume* volume 
  = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  //step size
  G4double stepSize = step->GetStepLength();  

  //collect energy step by step
  G4double edep = step->GetTotalEnergyDeposit();

  //track position coordinates
  G4StepPoint* pre = step->GetPreStepPoint();
  G4StepPoint* post = step->GetPostStepPoint();

  G4double xp_after,yp_after,zp_after,x_after,y_after,z_after,xp,yp,zp,x,y,z,xPix1,yPix1,zPix1,xPix,yPix,zPix,zOff;
  xp_after = pre->GetPosition().x();
  yp_after = pre->GetPosition().y();
  zp_after = pre->GetPosition().z();
  x_after = post->GetPosition().x();
  y_after = post->GetPosition().y();
  z_after = post->GetPosition().z();

  //geometry parameters
  G4double Det1SizeZ, Det2SizeZ, Dist, Strip1Depth, Strip1Length, Strip2Depth, Strip2Length, StripDist, StripWidth, posEndArm1, posBeginningArm2, BPIXSizeZ, pixelDepth, pixelX, pixelY, ROChor, ROCvert,  XangleDUT, XangleBPIX, YangleBPIX, ElField1, ElField2, zEndScint2, zlimMin, zlimMax;
  G4double posZAl, posZScint1, posZScint2, posZBPIX12, posZBPIX34, posZ1, posZ2, posZDUT;
  G4int div, div1, div2, divp1, divp2, dirx, diry, mod, ROC, row, col;
  G4int totalNbStrips, stripNo;
  Det1SizeZ=fDetectorconstruction->GetSize1();
  Det2SizeZ=fDetectorconstruction->GetSize2();
  Dist=fDetectorconstruction->GetDist();
  Strip1Depth=fDetectorconstruction->GetStrip1Depth();
  Strip1Length=fDetectorconstruction->GetStrip1Length();
  Strip2Depth=fDetectorconstruction->GetStrip2Depth();
  Strip2Length=fDetectorconstruction->GetStrip2Length();
  StripDist=fDetectorconstruction->GetStripDist();
  StripWidth=fDetectorconstruction->GetStripWidth();
  posEndArm1=-(fDetectorconstruction->Getpos_EndArm1Abs());
  posBeginningArm2=fDetectorconstruction->Getpos_BeginningArm2Abs();
  posZAl=fDetectorconstruction->GetPosZAl();
  posZScint1=fDetectorconstruction->GetPosZScint1();
  posZScint2=fDetectorconstruction->GetPosZScint2();
  posZBPIX12=fDetectorconstruction->GetPosZBPIX12();
  posZBPIX34=fDetectorconstruction->GetPosZBPIX34();
  posZ1=fDetectorconstruction->GetPosZ1();
  posZ2=fDetectorconstruction->GetPosZ2();
  posZDUT=fDetectorconstruction->GetPosZDUT();
  BPIXSizeZ=fDetectorconstruction->GetSizeBPIX();
  pixelX=fDetectorconstruction->GetPixelPitchX();
  pixelY=fDetectorconstruction->GetPixelPitchY();
  pixelDepth=fDetectorconstruction->GetPixelDepth();
  XangleDUT=fDetectorconstruction->GetDUTangleX();
  XangleBPIX=fDetectorconstruction->GetBPIXangleX();
  YangleBPIX=fDetectorconstruction->GetBPIXangleY();
  ElField1=fDetectorconstruction->GetElField1();
  ElField2=fDetectorconstruction->GetElField2();
  zEndScint2=fDetectorconstruction->Getpos_zEndScint2();
  zlimMin = zEndScint2 + 20.0*cm - stepSize/2;
  zlimMax = zEndScint2 + 20.0*cm + stepSize/2;
  totalNbStrips = 1016;
  ROChor = 52*pixelX + 2*pixelX;
  ROCvert = 80*pixelY + pixelY;

  //momentum direction when entering Pixel Detector 1 and exiting Pixel Detector 2
  G4double momDir1x, momDir1y, momDir1z, momDir2x, momDir2y, momDir2z;

  //position if the DUT wasn't rotated
  /*xp = xp_after;
  yp = yp_after*cos(XangleDUT) - zp_after*sin(XangleDUT);
  zp = yp_after*sin(XangleDUT) + zp_after*cos(XangleDUT);
  x = x_after;
  y = y_after*cos(XangleDUT) - z_after*sin(XangleDUT);
  z = y_after*sin(XangleDUT) + z_after*cos(XangleDUT);*/
  xp = xp_after;
  yp = yp_after;
  zp = zp_after;
  x = x_after;
  y = y_after;
  z = z_after;


  //position if the BPIX modules weren't rotated
  xPix1 = xp_after;
  yPix1 = yp_after;
  zPix1 = zp_after;
  xPix = xp_after;
  yPix = yp_after;
  zPix = zp_after;
  
  //secret number
  G4int iSecret, jSecret;  

  //get track length, track ID, track length, global time and track's vertex position of the current step
  G4Track* track = step->GetTrack();
  G4double length = track->GetTrackLength();
  G4double globalTime = track->GetGlobalTime();
  G4ThreeVector primVert = track->GetVertexPosition();
  G4int trackID = track->GetTrackID();

  //first, middle and last point of primary inside each detector
  G4ThreeVector pointAlent, pointAlex, pointScint1ent, pointScint1ex, pointScint2ent, pointScint2ex, pointPix1ent, pointPix1ex, pointPix2ent, pointPix2ex, pointDet1ent, pointDet1ex, pointDet2ent, pointDet2ex;
  G4ThreeVector pointPix1mid, pointPix2mid, pointDet1mid, pointDet2mid;

  //kinetic energy at entrance or exit of each detector
  G4double ekin, ekinPre, ekinPost, ekinAlent, ekinAlex, ekinScint1ent, ekinScint1ex, ekinScint2ent, ekinScint2ex, ekinPix1ent, ekinPix1ex, ekinPix2ent, ekinPix2ex, ekinDet1ent, ekinDet1ex, ekinDet2ent, ekinDet2ex;

  //point 20 cm after second scintillator
  G4ThreeVector point20Scint;

  //exit points after each volume
  G4ThreeVector pointExitAl, pointExitScint1, pointExitScint2, pointExitPix1, pointExitPix2, pointExitDet1, pointExitDet2;

  //particle definition
  const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();

  //track length of primary particle
  if (track->GetTrackID() == 1)  {
    fRunaction->AddTrackLength(stepSize);
    G4AnalysisManager::Instance()->FillH1(17,stepSize);

    fEventaction->AddPrimaryTrackLength(stepSize);

    ekinPre = pre->GetKineticEnergy(); //post?
    ekinPost = pre->GetKineticEnergy(); //post?
    /*G4cout
       << "\n ekinPre = "
       << G4BestUnit(ekinPre, "Energy")
       << "\n ekinPost = "
       << G4BestUnit(ekinPost, "Energy")
       << "\n z = "
       << G4BestUnit(z, "Length")
       << G4endl;*/
    zOff = z + 23.5*cm;
    G4AnalysisManager::Instance()->FillH2(10,zOff,ekinPre);

    
    //kinetic energy of primary particle at entrance and exit of each volume
    if (volume == fDetectorconstruction->GetAl())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	   ekinAlent = pre->GetKineticEnergy();
        }

	if (post->GetStepStatus() == fGeomBoundary)   {
	   ekinAlex = post->GetKineticEnergy();
	   fEventaction->AddEnergyExAl(ekinAlex);
	   fEventaction->AddPointExitAl(post->GetPosition());

        }
    }

    if (volume == fDetectorconstruction->GetScint1())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	   ekinScint1ent = pre->GetKineticEnergy();
	   fEventaction->AddEnergyEntScintillator1(ekinScint1ent);
	   fEventaction->AddTimeEntranceScint1(globalTime);
        }

	if (post->GetStepStatus() == fGeomBoundary)   {
	   ekinScint1ex = post->GetKineticEnergy();
	   fEventaction->AddEnergyExScintillator1(ekinScint1ex);
	   fEventaction->AddPointExitScint1(post->GetPosition());

        }
    }

    if (volume == fDetectorconstruction->GetBPIX12Sen())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    pointPix1ent = pre->GetPosition();
	    momDir1x = track->GetMomentumDirection().x();
	    momDir1y = track->GetMomentumDirection().y();
	    momDir1z = track->GetMomentumDirection().z();

            fEventaction->AddMomentumDirection1(momDir1x, momDir1y, momDir1z);

    	    fEventaction->AddPointPix1entReal(pointPix1ent);

            //G4cout 
               //<< "\n Primary has entered BPIX 12 at: " 
               //<< G4BestUnit(pointPix1ent, "Length")
               //<< "\n with momentum direction: " 
               //<< momDir1x
               //<< " "
               //<< momDir1y
               //<< " "
               //<< momDir1z
               //<< G4endl;

	   ekinPix1ent = pre->GetKineticEnergy();
	   fEventaction->AddEnergyEntPixelLayer1(ekinPix1ent);
        }
    }
    if (volume == fDetectorconstruction->GetBPIX12())   {
	if (post->GetStepStatus() == fGeomBoundary)   {
	    pointPix1ex = post->GetPosition();

            //G4cout 
               //<< "\n Primary has exited Pixel Detector 1 at: " 
               //<< G4BestUnit(pointPix1ex, "Length")
               //<< G4endl;

	   ekinPix1ex = post->GetKineticEnergy();
	   fEventaction->AddEnergyExPixelLayer1(ekinPix1ex);
	   fEventaction->AddPointExitPixel12(post->GetPosition());
        }
    }

    if (volume == fDetectorconstruction->GetDet1())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    pointDet1ent = pre->GetPosition();

    	    fEventaction->AddPointDet1entReal(pointDet1ent);

            //G4cout 
               //<< "\n Primary has entered Strip Detector 1 at: " 
               //<< G4BestUnit(pointDet1ent, "Length")
               //<< G4endl;

	   ekinDet1ent = pre->GetKineticEnergy();
	   fEventaction->AddEnergyEntSensor1(ekinDet1ent);
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    pointDet1ex = post->GetPosition();

            //G4cout 
               //<< "\n Primary has exited Strip Detector 1 at: " 
               //<< G4BestUnit(pointDet1ex, "Length")
               //<< G4endl;

	   ekinDet1ex = post->GetKineticEnergy();
	   fEventaction->AddEnergyExSensor1(ekinDet1ex);
	   fEventaction->AddPointExitSensor1(post->GetPosition());
        }
    }

    if (volume == fDetectorconstruction->GetDet2())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    pointDet2ent = pre->GetPosition();

    	    fEventaction->AddPointDet2entReal(pointDet2ent);

            //G4cout 
               //<< "\n Primary has entered Strip Detector 2 at: " 
               //<< G4BestUnit(pointDet2ent, "Length")
               //<< G4endl;

	   ekinDet2ent = pre->GetKineticEnergy();
	   fEventaction->AddEnergyEntSensor2(ekinDet2ent);
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    pointDet2ex = post->GetPosition();

            //G4cout 
               //<< "\n Primary has exited Strip Detector 2 at: " 
               //<< G4BestUnit(pointDet2ex, "Length")
               //<< G4endl;

	   ekinDet2ex = post->GetKineticEnergy();
	   fEventaction->AddEnergyExSensor2(ekinDet2ex);
	   fEventaction->AddPointExitSensor2(post->GetPosition());
        }
    }

    if (volume == fDetectorconstruction->GetBPIX34Sen())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    pointPix2ent = pre->GetPosition();

    	    fEventaction->AddPointPix2entReal(pointPix2ent);

            //G4cout 
               //<< "\n Primary has entered Pixel Detector 2 at: " 
               //<< G4BestUnit(pointPix2ent, "Length")
               //<< G4endl;

	   ekinPix2ent = pre->GetKineticEnergy();
	   fEventaction->AddEnergyEntPixelLayer2(ekinPix2ent);
        }
    }
    if (volume == fDetectorconstruction->GetBPIX34())   {
	if (post->GetStepStatus() == fGeomBoundary)   {
	    pointPix2ex = post->GetPosition();
	    momDir2x = track->GetMomentumDirection().x();
	    momDir2y = track->GetMomentumDirection().y();
	    momDir2z = track->GetMomentumDirection().z();

            fEventaction->AddMomentumDirection2(momDir2x, momDir2y, momDir2z);

            //G4cout 
               //<< "\n Primary has exited Pixel Detector 2 at: " 
               //<< G4BestUnit(pointPix2ex, "Length")
               //<< "\n with momentum direction: " 
               //<< momDir2x
               //<< " "
               //<< momDir2y
               //<< " "
               //<< momDir2z
               //<< G4endl;

	   ekinPix2ex = post->GetKineticEnergy();
	   fEventaction->AddEnergyExPixelLayer2(ekinPix2ex);
	   fEventaction->AddPointExitPixel34(post->GetPosition());
        }
    }
    /*G4cout 
           << "\n zp_after = " 
           << G4BestUnit(zp_after, "Length")
           << "\n zlimMin = " 
           << G4BestUnit(zlimMin, "Length")
           << "\n zlimMax = " 
           << G4BestUnit(zlimMax, "Length")
           << G4endl;*/

    if (volume == fDetectorconstruction->GetScint2())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	   ekinScint2ent = pre->GetKineticEnergy();
	   fEventaction->AddEnergyEntScintillator2(ekinScint2ent);
	   fEventaction->AddTimeEntranceScint2(globalTime);
        }

	if (post->GetStepStatus() == fGeomBoundary)   {
	   ekinScint2ex = post->GetKineticEnergy();
	   fEventaction->AddEnergyExScintillator2(ekinScint2ex);
	   fEventaction->AddPointExitScint2(post->GetPosition());
        }
    }

    if ((zp_after <= zlimMax) && (zp_after >= zlimMin))  {
    //zlimMin = zEndScint2 + 20.0*cm;
    /*G4cout 
           << "\n zlimMin = " 
           << G4BestUnit(zlimMin, "Length")
           << G4endl;*/
    //if (zp_after >= zEndScint2)  {
        point20Scint = pre->GetPosition();
        fEventaction->AddPointAfterScint2(point20Scint);
        /*G4cout 
           << "\n zp_after = " 
           << G4BestUnit(zp_after, "Length")
           << "\n zEndScint2 = " 
           << G4BestUnit(zEndScint2, "Length")
           << "\n point20Scint = " 
           << G4BestUnit(point20Scint, "Length")
           << G4endl;*/

    }
  }

  //track length of secondaries and tertiaries calculation
  if (track->GetParentID() == 1)  {
    fRunaction->AddSecTrackLength(stepSize);
    G4AnalysisManager::Instance()->FillH1(18,stepSize);

    fEventaction->AddSecondaryTrackLength(stepSize);

    if ((track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet1()) || (track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet2()))  {
         fRunaction->AddSecDetTrackLength(stepSize);
         G4AnalysisManager::Instance()->FillH1(20,stepSize);

         fEventaction->AddSecondaryDetTrackLength(stepSize);
         fEventaction->TrackCheck(length);

         //G4cout 
            //<< "\n Secondary produced at: " 
            //<< G4BestUnit(primVert, "Length")
            //<< "\n with track ID: " 
            //<< trackID
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector has a track length: " 
            //<< G4BestUnit(length, "Length")
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector's current position: " 
            //<< G4BestUnit(x, "Length")
            //<< G4BestUnit(y, "Length")
            //<< G4BestUnit(z, "Length")
            //<< G4endl;
    }

    if (track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet1())  {
         fRunaction->AddSecDetTrackLength(stepSize);
         G4AnalysisManager::Instance()->FillH1(20,stepSize);

         fEventaction->AddSecondaryDetTrackLength(stepSize);
         fEventaction->TrackCheck(length);

         //G4cout 
            //<< "\n Secondary produced at: " 
            //<< G4BestUnit(primVert, "Length")
            //<< "\n with track ID: " 
            //<< trackID
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside Detector 1 has a track length: " 
            //<< G4BestUnit(length, "Length")
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector's current position: " 
            //<< G4BestUnit(x, "Length")
            //<< G4BestUnit(y, "Length")
            //<< G4BestUnit(z, "Length")
            //<< G4endl;
    }

    if (track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet2())  {
         fRunaction->AddSecDetTrackLength(stepSize);
         G4AnalysisManager::Instance()->FillH1(20,stepSize);

         fEventaction->AddSecondaryDetTrackLength(stepSize);
         fEventaction->TrackCheck(length);

         //G4cout 
            //<< "\n Secondary produced at: " 
            //<< G4BestUnit(primVert, "Length")
            //<< "\n with track ID: " 
            //<< trackID
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside Detector 2 has a track length: " 
            //<< G4BestUnit(length, "Length")
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector's current position: " 
            //<< G4BestUnit(x, "Length")
            //<< G4BestUnit(y, "Length")
            //<< G4BestUnit(z, "Length")
            //<< G4endl;
    }

    if ((track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet1()) && (volume == fDetectorconstruction->GetDet2()))  {
         //G4cout 
            //<< "\n Secondary produced inside Detector 1 has reached Detector 2. TrackID = " 
            //<< trackID
            //<< G4endl;

         if(particleDefinition == G4Positron::Definition())   {
              //G4cout 
                 //<< "\n It was a positron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Electron::Definition())   {
              //G4cout 
                 //<< "\n It was an electron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Gamma::Definition())   {
              //G4cout 
                 //<< "\n It was a photon" 
                 //<< G4endl;
         }


   	 if ((z >= posZDUT + Dist/2 + Det2SizeZ - Strip2Depth) && (z <= (posZDUT + Dist/2 + Det2SizeZ))&& (x >= (-Strip2Length/2)) && (x <= Strip2Length/2)) {
     	  div = y/(StripWidth+StripDist);
	  if (y == 0)   {
	     iSecret = rand() % 99;
	     if (iSecret < 50)   {
		  stripNo = totalNbStrips/2 + div;
	     }
	     if (iSecret >= 50)   {
	  	  stripNo = totalNbStrips/2 + 1 + div;
	     }
	  }
     	  if (y > 0)    {
	     stripNo = totalNbStrips/2 + 1 + div;
     	  }
     	  if (y < 0)    {
             stripNo = totalNbStrips/2 + div;
     	  }
     	  if ((stripNo <= totalNbStrips) && (stripNo > 0))   {
     	     //G4cout
	         //<< "\n The secondary has passed below strip No "
	         //<< stripNo
                 //<< G4endl;
          }
         }
    }

    if ((track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet2()) && (volume == fDetectorconstruction->GetDet1()))  {
         //G4cout 
            //<< "\n Secondary produced inside Detector 2 has reached Detector 1. TrackID = " 
            //<< trackID
            //<< G4endl;

         if(particleDefinition == G4Positron::Definition())   {
              //G4cout 
                 //<< "\n It was a positron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Electron::Definition())   {
              //G4cout 
                 //<< "\n It was an electron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Gamma::Definition())   {
              //G4cout 
                 //<< "\n It was a photon" 
                 //<< G4endl;
         }

    	 if ((z >= posZDUT - Dist/2 - Det1SizeZ) && (z <= posZDUT - Dist/2 - Det1SizeZ + Strip1Depth) && (x >= (-Strip1Length/2)) && (x <= Strip1Length/2)) {
     	  div = y/(StripWidth+StripDist);
	  if (y == 0)   {
	     iSecret = rand() % 99;
	     if (iSecret < 50)   {
		  stripNo = totalNbStrips/2 + div;
	     }
	     if (iSecret >= 50)   {
	       	  stripNo = totalNbStrips/2 + 1 + div;
	     }
	  }
     	  if (y > 0)    {
	      stripNo = totalNbStrips/2 + 1 + div;
     	  }
     	  if (y < 0)    {
              stripNo = totalNbStrips/2 + div;
     	  }
     	  if ((stripNo <= totalNbStrips) && (stripNo > 0))  {
     	      //G4cout
	   	  //<< "\n The secondary has passed below strip No "
	   	  //<< stripNo
           	  //<< G4endl;
     	  }
    	 }
    }
  }

  if (track->GetParentID() > 1)  {
    fRunaction->AddTertTrackLength(stepSize);
    G4AnalysisManager::Instance()->FillH1(19,stepSize);

    fEventaction->AddTertiaryTrackLength(stepSize);
  }

 //continuous energy deposit per event  

 if (volume == fDetectorconstruction->GetBPIX12Sen()) {

   if ((zPix >= (posZBPIX12 - BPIXSizeZ/2))  &&  (zPix < (posZBPIX12 - BPIXSizeZ/2 + pixelDepth))) {
     fEventaction->AddEnergyDepositL1(edep);

     mod = 1; //module ID
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp2 = yPix/pixelY;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if ((diry == 0) && (dirx == 0))   {
         if (div1 == 3)  {
             ROC = 1;
         }
         if (div1 == 2)  {
             ROC = 2;
         }
         if (div1 == 1)  {
             ROC = 3;
         }
         if (div1 == 0)  {
             ROC = 4;
         }
     }
     if ((diry == 0) && (dirx == 1))   {
         if (div1 == 0)  {
             ROC = 5;
         }
         if (div1 == -1)  {
             ROC = 6;
         }
         if (div1 == -2)  {
             ROC = 7;
         }
         if (div1 == -3)  {
             ROC = 8;
         }
     }
     if ((diry == 1) && (dirx == 0))   {
         if (div1 == 3)  {
             ROC = 9;
         }
         if (div1 == 2)  {
             ROC = 10;
         }
         if (div1 == 1)  {
             ROC = 11;
         }
         if (div1 == 0)  {
             ROC = 12;
         }
     }
     if ((diry == 1) && (dirx == 1))   {
         if (div1 == 0)  {
             ROC = 13;
         }
         if (div1 == -1)  {
             ROC = 14;
         }
         if (div1 == -2)  {
             ROC = 15;
         }
         if (div1 == -3)  {
             ROC = 16;
         }
     }
     if (diry == 0)  {
	 if ((divp2 == 0) || (divp2 == 1))  {
	    row = 80;
         }
         if ((divp2 > 1) && (divp2 <= 80))  {
	    row = 80 - divp2 + 1;
         } 
     }
     if (diry == 1)  {
	 if ((divp2 == -79) || (divp2 == -80))  {
	    row = 80;
         }
         if ((divp2 > -79) && (divp2 <= 0))  {
	    row = -divp2 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp1 = (xPix-3*ROChor)/pixelX;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp1 = (xPix-2*ROChor)/pixelX;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp1 = (xPix-1*ROChor)/pixelX;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp1 = (xPix-0*ROChor)/pixelX;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp1 = (xPix+0*ROChor)/pixelX;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp1 = (xPix+1*ROChor)/pixelX;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp1 = (xPix+2*ROChor)/pixelX;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp1 = (xPix+3*ROChor)/pixelX;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp1 == 0) || (divp1 == 1))  {
	   col = 52;
	}
	if ((divp1 == 52) || (divp1 == 53))  {
	   col = 1;
	}
	if ((divp1 < 52) && (divp1 > 1))  {
	   col = 52 - divp1 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp1 == 0) || (divp1 == -1))  {
	   col = 1;
	}
	if ((divp1 == -52) || (divp1 == -53))  {
	   col = 52;
	}
	if ((divp1 > -52) && (divp1 < -1))  {
	   col = -divp1;
	}
     }

     if ((mod>=1) && (mod<3) && (ROC>=1) && (ROC<17) && (row>=1) && (row<81) && (col>=1) && (col<53))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
     }
   }
 }
 if (volume == fDetectorconstruction->GetBPIX12()) {
   if (zPix >= (posZBPIX12 - BPIXSizeZ/2 + pixelDepth)) {
        fEventaction->AddEnergyDepositROCL1(edep);
   }
 }

 if (volume == fDetectorconstruction->GetBPIX34Sen()) {
   if ((zPix >= posZBPIX34 - BPIXSizeZ/2)  &&  (zPix < (posZBPIX34 - BPIXSizeZ/2 + pixelDepth))) {
     fEventaction->AddEnergyDepositL2(edep);

     mod = 2; //module ID
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp2 = yPix/pixelY;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if ((diry == 0) && (dirx == 0))   {
         if (div1 == 3)  {
             ROC = 1;
         }
         if (div1 == 2)  {
             ROC = 2;
         }
         if (div1 == 1)  {
             ROC = 3;
         }
         if (div1 == 0)  {
             ROC = 4;
         }
     }
     if ((diry == 0) && (dirx == 1))   {
         if (div1 == 0)  {
             ROC = 5;
         }
         if (div1 == -1)  {
             ROC = 6;
         }
         if (div1 == -2)  {
             ROC = 7;
         }
         if (div1 == -3)  {
             ROC = 8;
         }
     }
     if ((diry == 1) && (dirx == 0))   {
         if (div1 == 3)  {
             ROC = 9;
         }
         if (div1 == 2)  {
             ROC = 10;
         }
         if (div1 == 1)  {
             ROC = 11;
         }
         if (div1 == 0)  {
             ROC = 12;
         }
     }
     if ((diry == 1) && (dirx == 1))   {
         if (div1 == 0)  {
             ROC = 13;
         }
         if (div1 == -1)  {
             ROC = 14;
         }
         if (div1 == -2)  {
             ROC = 15;
         }
         if (div1 == -3)  {
             ROC = 16;
         }
     }
     if (diry == 0)  {
	 if ((divp2 == 0) || (divp2 == 1))  {
	    row = 80;
         }
         if ((divp2 > 1) && (divp2 <= 80))  {
	    row = 80 - divp2 + 1;
         } 
     }
     if (diry == 1)  {
	 if ((divp2 == -79) || (divp2 == -80))  {
	    row = 80;
         }
         if ((divp2 > -79) && (divp2 <= 0))  {
	    row = -divp2 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp1 = (xPix-3*ROChor)/pixelX;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp1 = (xPix-2*ROChor)/pixelX;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp1 = (xPix-1*ROChor)/pixelX;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp1 = (xPix-0*ROChor)/pixelX;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp1 = (xPix+0*ROChor)/pixelX;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp1 = (xPix+1*ROChor)/pixelX;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp1 = (xPix+2*ROChor)/pixelX;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp1 = (xPix+3*ROChor)/pixelX;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp1 == 0) || (divp1 == 1))  {
	   col = 52;
	}
	if ((divp1 == 52) || (divp1 == 53))  {
	   col = 1;
	}
	if ((divp1 < 52) && (divp1 > 1))  {
	   col = 52 - divp1 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp1 == 0) || (divp1 == -1))  {
	   col = 1;
	}
	if ((divp1 == -52) || (divp1 == -53))  {
	   col = 52;
	}
	if ((divp1 > -52) && (divp1 < -1))  {
	   col = -divp1;
	}
     } 

     if ((mod>=1) && (mod<3) && (ROC>=1) && (ROC<17) && (row>=1) && (row<81) && (col>=1) && (col<53))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
     }
   }
 }
 if (volume == fDetectorconstruction->GetBPIX34()) {
   if (zPix >= (posZBPIX34 - BPIXSizeZ/2 + pixelDepth)) {
        fEventaction->AddEnergyDepositROCL2(edep);
   }
 }


 if (volume == fDetectorconstruction->GetDet1()) {
   fEventaction->AddEnergyDeposit1(edep);

   if ((z >= posZDUT - Dist/2 - Det1SizeZ) && (z <= posZDUT - Dist/2 - Det1SizeZ + Strip1Depth) && (x >= (-Strip1Length/2)) && (x <= Strip1Length/2)) {
     div = y/(StripWidth+StripDist);
     if (y == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   stripNo = totalNbStrips/2 + div;
     	}
     	if (iSecret >= 50)   {
   	   stripNo = totalNbStrips/2 + 1 + div;
     	}
     }
     if (y > 0)    {
	stripNo = totalNbStrips/2 + 1 + div;
     }
     if (y < 0)    {
        stripNo = totalNbStrips/2 + div;
     }
     if ((stripNo <= totalNbStrips) && (stripNo > 0))  {
     	/*G4cout
	   << "\n Continuous energy deposition below strip No "
           << stripNo
           << "\n y = "
	   << y
           << "\n div = "
	   << div
           << G4endl;*/
        if (x > 0)  {
	   fEventaction->AddEnergyStrip1a(edep,stripNo);
        }
        if (x < 0)  {
	   fEventaction->AddEnergyStrip1b(edep,stripNo);
        }
        if (x == 0)  {
	    jSecret = rand() % 99;
	    if (jSecret < 50)   {
	       fEventaction->AddEnergyStrip1a(edep,stripNo);
	    }
	    if (jSecret >= 50)   {
	       fEventaction->AddEnergyStrip1b(edep,stripNo);
	    }
        }
     }
   }
 }

 if (volume == fDetectorconstruction->GetDet2()) {
   fEventaction->AddEnergyDeposit2(edep);

   if ((z >= posZDUT + Dist/2 + Det2SizeZ - Strip2Depth) && (z <= (posZDUT + Dist/2 + Det2SizeZ))&& (x >= (-Strip2Length/2)) && (x <= Strip2Length/2)) {
     div = y/(StripWidth+StripDist);
     if (y == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   stripNo = totalNbStrips/2 + div;
     	}
     	if (iSecret >= 50)   {
   	   stripNo = totalNbStrips/2 + 1 + div;
     	}
     }
     if (y > 0)    {
	stripNo = totalNbStrips/2 + 1 + div;
     }
     if (y < 0)    {
        stripNo = totalNbStrips/2 + div;
     }
     if ((stripNo <= totalNbStrips)&& (stripNo > 0))   {
     	//G4cout
	   //<< "\n Continuous energy deposition below strip No "
           //<< stripNo
           //<< "\n y = "
	   //<< y
           //<< "\n div = "
	   //<< div
           //<< G4endl;
        if (x > 0)  {
	   fEventaction->AddEnergyStrip2a(edep,stripNo);
        }
        if (x < 0)  {
	   fEventaction->AddEnergyStrip2b(edep,stripNo);
        }
        if (x == 0)  {
	    jSecret = rand() % 99;
	    if (jSecret < 50)   {
	       fEventaction->AddEnergyStrip2a(edep,stripNo);
	    }
	    if (jSecret >= 50)   {
	       fEventaction->AddEnergyStrip2b(edep,stripNo);
	    }
        }
     }
   }
 }

 if (volume == fDetectorconstruction->GetScint1()) {
   fEventaction->AddEnergyDepositSc1(edep);
 }

 if (volume == fDetectorconstruction->GetScint2()) {
   fEventaction->AddEnergyDepositSc2(edep);
 }

 if (volume == fDetectorconstruction->GetAl()) {
   fEventaction->AddEnergyDepositAl(edep);
 }

 if ((volume == fDetectorconstruction->GetWorld()) && (volume != fDetectorconstruction->GetDet1()) && (volume != fDetectorconstruction->GetDet2()) && (volume != fDetectorconstruction->GetBPIX34()) && (volume != fDetectorconstruction->GetScint1()) && (volume != fDetectorconstruction->GetScint2()) && (volume != fDetectorconstruction->GetAl())){
   fEventaction->AddEnergyDepositAir(edep); // || (volume == fDetectorconstruction->GetBPIX12()))
 }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

