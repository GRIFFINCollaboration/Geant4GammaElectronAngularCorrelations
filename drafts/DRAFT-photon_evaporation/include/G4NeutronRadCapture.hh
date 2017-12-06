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
// $Id: G4NeutronRadCapture.hh 90590 2015-06-04 13:37:40Z gcosmo $
//
//
// Geant4 header : G4NeutronRadCapture
// Created:  31 August 2009
// Author  V.Ivanchenko
//  
// Modified:
//
// Class Description
// Sampling of neutron radiative capture 
// Class Description - End
//

#ifndef G4NeutronRadCapture_h
#define G4NeutronRadCapture_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"

//const G4ParticleDefinition* theDef;

class G4PhotonEvaporation;
class G4IonTable;

class G4NeutronRadCapture : public G4HadronicInteraction
{
public:

  G4NeutronRadCapture();

  virtual ~G4NeutronRadCapture();
 
  G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
				  G4Nucleus & targetNucleus);

private:

  G4NeutronRadCapture & operator=(const G4NeutronRadCapture &right);
  G4NeutronRadCapture(const G4NeutronRadCapture&);

  G4double lowestEnergyLimit;
  G4double minExcitation;
  G4PhotonEvaporation* photonEvaporation;
  G4IonTable*  theTableOfIons;
  G4LorentzVector lab4mom;

};

#endif
