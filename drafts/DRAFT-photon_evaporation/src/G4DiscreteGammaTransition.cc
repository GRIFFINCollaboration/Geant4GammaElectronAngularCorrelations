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
// $Id: G4DiscreteGammaTransition.cc 87257 2014-11-28 08:02:05Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4DiscreteGammaTransition
//
//      Author:        Maria Grazia Pia   (pia@genova.infn.it)
//
//      Creation date: 23 October 1998
//
//      Modifications:
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              i) added G4int _nucleusZ initialise it through the constructor
//              ii) modified SelectGamma() to allow the generation of
//                  conversion e-
//              iii) added #include G4AtomicShells.hh
//
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added renormalization to determine whether transition leads to
//              electron or gamma in SelectGamma()
//
//        19 April 2010, J. M. Quesada.
//              Corrections added for taking into account mismatch between
//              tabulated gamma energies and level energy differences
//              (fake photons eliminated)
//
//        9 May 2010, V.Ivanchenko
//              Removed unphysical corretions of gamma energy; fixed default
//              particle as gamma; do not subtract bounding energy in case of
//              electron emmision
//
//	  3 November 2011, L. Desorgher
//		Extend the use of the code for Z>100 by not calling
//              G4AtomicShells::GetBindingEnergy for Z>100
//		For Z>100 the binding energy is set to 0 and atomic
//              relaxation is not simulated in G4RadDecay
//
// -------------------------------------------------------------------

#include "G4DiscreteGammaTransition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4AtomicShells.hh"
#include "G4NuclearLevel.hh"
#include "G4NuclearLevelStore.hh"
#include "G4Pow.hh"
#include "G4Log.hh"

static const G4double tolerance = 10*CLHEP::keV;

G4DiscreteGammaTransition::G4DiscreteGammaTransition(
  const G4NuclearLevel* level, G4int Z, G4int verb)
  : orbitE(-1), bondE(0.),  gammaEnergy(0.),  excitation(0.),
    gammaCreationTime(0.), aGamma(true), icm(false)
{
  verbose = verb;
  Update(level, Z);
}

G4DiscreteGammaTransition::~G4DiscreteGammaTransition()
{}

void G4DiscreteGammaTransition::SelectGamma()
{
  // default gamma
  aGamma = true;
  gammaEnergy = 0.;

  G4int nGammas = aLevel->NumberOfGammas();
  if (nGammas > 0) {
    G4int iGamma = 0;
    if(1 < nGammas) {
      G4double random = G4UniformRand();

      //G4cout << "G4DiscreteGammaTransition::SelectGamma  N= "
      //	     << nGammas << " rand= " << random << G4endl;
      for(iGamma=0; iGamma<nGammas; ++iGamma) {
    //G4cout << iGamma << "  prob= "
    //   << (aLevel->GammaCumulativeProbabilities())[iGamma] << G4endl;
    if(random <= (aLevel->GammaCumulativeProbabilities())[iGamma])
      { break; }
      }
    }
    /*
    G4cout << "Elevel(MeV)= " << aLevel->Energy()
       << " Etran(MeV)= " << (aLevel->GammaEnergies())[iGamma]
       << " Eexc(MeV)= " << excitation << G4endl;
    */
    // VI 2014: initial excitation energy may be not exactly energy of the level
    //          final excitation energy is always energy of some level or zero
    //          transition to the ground state should be always equal to
    //          the excitation energy
    gammaEnergy = (aLevel->GammaEnergies())[iGamma]
      + excitation - aLevel->Energy();

    // this check is needed to remove cases when nucleaus is left in
    // slightly excited state which will require very low energy
    // gamma emission
    if(excitation <= gammaEnergy + tolerance) { gammaEnergy = excitation; }
    //JMQ:
    //1)If chosen gamma energy is close enough to excitation energy,
    //  the later is used instead for gamma dacey to gs (it guarantees
    //  energy conservation)
    //2)For energy conservation, level energy differences instead of
    //  tabulated gamma energies must be used (origin of final fake photons)

    // VI: remove fake photons - applied only for the last transition
    //     do not applied on each transition
    //if(std::fabs(excitation - gammaEnergy) < tolerance) {
    //  gammaEnergy = excitation;
    //}

    //  JMQ: Warning: the following check is needed to avoid loops:
    //  Due essentially to missing nuclear levels in data files, it is
    //  possible that gammaEnergy is so low as the nucleus doesn't change
    //  its level after the transition.
    //  When such case is found, force the full deexcitation of the nucleus.
    //
    //    NOTE: you should force the transition to the next lower level,
    //          but this change needs a more complex revision of actual
    //          design.
    //          I leave this for a later revision.

    // VI: the check is needed to remove very low-energy gamma
    if (gammaEnergy < tolerance) { gammaEnergy = excitation; }
    /*
    G4cout << "G4DiscreteGammaTransition::SelectGamma: " << gammaEnergy
       << " Eexc= " << excitation
       << " icm: " << icm << G4endl;
    */

    // Rishita -- TODO: add vectors for other three cases
    // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    _levelEnergy         = (aLevel->Energy());
    
    _a2_gg               = (aLevel->A2_GG())[iGamma];
    _a4_gg               = (aLevel->A4_GG())[iGamma];
    _a6_gg               = (aLevel->A6_GG())[iGamma];
    _a8_gg               = (aLevel->A8_GG())[iGamma];
    _a10_gg              = (aLevel->A10_GG())[iGamma];
    
    _a2_ge               = (aLevel->A2_GE())[iGamma];
    _a4_ge               = (aLevel->A4_GE())[iGamma];
    _a6_ge               = (aLevel->A6_GE())[iGamma];
    _a8_ge               = (aLevel->A8_GE())[iGamma];
    _a10_ge              = (aLevel->A10_GE())[iGamma];
    
    _a2_eg               = (aLevel->A2_EG())[iGamma];
    _a4_eg               = (aLevel->A4_EG())[iGamma];
    _a6_eg               = (aLevel->A6_EG())[iGamma];
    _a8_eg               = (aLevel->A8_EG())[iGamma];
    _a10_eg              = (aLevel->A10_EG())[iGamma];
    
    _a2_ee               = (aLevel->A2_EE())[iGamma];
    _a4_ee               = (aLevel->A4_EE())[iGamma];
    _a6_ee               = (aLevel->A6_EE())[iGamma];
    _a8_ee               = (aLevel->A8_EE())[iGamma];
    _a10_ee              = (aLevel->A10_EE())[iGamma];
    
    _maxWTheta_gg           = (aLevel->MaxWTheta_GG())[iGamma];
    _maxWTheta_ge           = (aLevel->MaxWTheta_GE())[iGamma];
    _maxWTheta_eg           = (aLevel->MaxWTheta_EG())[iGamma];
    _maxWTheta_ee           = (aLevel->MaxWTheta_EE())[iGamma];
    
    _higherLevelEnergy   = (aLevel->HigherLevelEnergy())[iGamma];
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    // now decide whether Internal Coversion electron should be emitted instead
    if (icm) {
      G4double random = G4UniformRand();
      if ( random <= (aLevel->TotalConvertionProbabilities())[iGamma]
       *(aLevel->GammaWeights())[iGamma]
       /((aLevel->TotalConvertionProbabilities())[iGamma]
         *(aLevel->GammaWeights())[iGamma]
         +(aLevel->GammaWeights())[iGamma]))
    {
      G4int iShell = 9;
      random = G4UniformRand() ;
      if ( random <= (aLevel->KConvertionProbabilities())[iGamma])
        { iShell = 0;}
      else if ( random <= (aLevel->L1ConvertionProbabilities())[iGamma])
        { iShell = 1;}
      else if ( random <= (aLevel->L2ConvertionProbabilities())[iGamma])
        { iShell = 2;}
      else if ( random <= (aLevel->L3ConvertionProbabilities())[iGamma])
        { iShell = 3;}
      else if ( random <= (aLevel->M1ConvertionProbabilities())[iGamma])
        { iShell = 4;}
      else if ( random <= (aLevel->M2ConvertionProbabilities())[iGamma])
        { iShell = 5;}
      else if ( random <= (aLevel->M3ConvertionProbabilities())[iGamma])
        { iShell = 6;}
      else if ( random <= (aLevel->M4ConvertionProbabilities())[iGamma])
        { iShell = 7;}
      else if ( random <= (aLevel->M5ConvertionProbabilities())[iGamma])
        { iShell = 8;}
      // the following is needed to match the ishell to that used in
      // G4AtomicShells
      /*
      if ( iShell == 9) {
        if ( (nucleusZ < 28) && (nucleusZ > 20)) {
          iShell--;
        } else if ( nucleusZ == 20 || nucleusZ == 19 ) {
          iShell = iShell -2;
        }
      }
      */
      //L.Desorgher 02/11/2011
      //Atomic shell information is available in Geant4 only up top Z=100
      //To extend the photo evaporation code to Z>100  the call
      // to G4AtomicShells::GetBindingEnergy should be forbidden for Z>100
      bondE = 0.;
      if (nucleusZ <=100) {
        bondE = G4AtomicShells::GetBindingEnergy(nucleusZ, iShell);
      }
      if (verbose > 1) {
        G4cout << "G4DiscreteGammaTransition: nucleusZ = " <<nucleusZ
           << " , iShell = " << iShell
           << " , Shell binding energy = " << bondE/keV
           << " keV " << G4endl;
      }

      // last check on energy
      if(gammaEnergy >  bondE + tolerance) {
        orbitE = iShell;
        aGamma = false ;   // emitted is not a gamma now
    //    gammaEnergy -= bondE;
      }
      //G4cout << "gammaEnergy = " << gammaEnergy << G4endl;
    }
    }

    G4double tau = aLevel->HalfLife() / G4Pow::GetInstance()->logZ(2);

    //09.05.2010 VI rewrite samling of decay time
    //              assuming ordinary exponential low
    gammaCreationTime = 0.;
    if(tau > 0.0) {  gammaCreationTime = -tau*G4Log(G4UniformRand()); }
  }
  //G4cout << "G4DiscreteGammaTransition end nGamma= " << nGammas
  //	 << "  Egamma= " << gammaEnergy << G4endl;
}

G4double G4DiscreteGammaTransition::GetGammaEnergy()
{
  return gammaEnergy;
}

G4double G4DiscreteGammaTransition::GetGammaCreationTime()
{
  return gammaCreationTime;
}

void G4DiscreteGammaTransition::SetEnergyFrom(G4double energy)
{
  excitation = energy;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G4double G4DiscreteGammaTransition::GetThetaFromWTheta(G4double higherLevelEnergy, G4int casenumber)
{
  G4double  value, WTheta_gg, WTheta_ge, WTheta_eg, WTheta_ee, theta, normy_gg, normy_ge, normy_eg, normy_ee;
  G4bool    boolhlevelindex, foundTheta;
  G4int     hlevelindex;

//G4cout << "------------ G4DiscreteGammaTransition START -----------------------------------------------" << G4endl;

  if(higherLevelEnergy<0) {
    value = std::acos(1. - 2. * G4UniformRand());
  }
  else {
    boolhlevelindex = false;
    for (G4int i=0; i<_higherLevelEnergy.size(); i++) {
      if(_higherLevelEnergy[i] == higherLevelEnergy) {
        boolhlevelindex = true;
        hlevelindex = i;
        break;
      }
    }
    if(!boolhlevelindex) {
      G4cout << "higherLevelEnergy = " << higherLevelEnergy/keV << " keV" << G4endl;
      G4cout << "Didn't find this higher level index" << G4endl;
      exit(1);
    }
// Rishita
/*    std::vector<const G4ParticleDefinition*>* rep = G4NuclearLevelStore::GetInstance()->GetResults();
    if (!rep) {G4cout << "none" <<G4endl; } 
    else{for (G4int i = 0; i < rep->size(); i++) {
        const G4ParticleDefinition* n = rep->at(i);
        G4String name = n->GetParticleName();
        G4cout << "Name :" << name << G4endl;
    }}
*/
    foundTheta = false;
    while (!foundTheta) {
      theta = M_PI*G4UniformRand();

      normy_gg = G4UniformRand()*_maxWTheta_gg[hlevelindex];
      normy_ge = G4UniformRand()*_maxWTheta_ge[hlevelindex];
      normy_eg = G4UniformRand()*_maxWTheta_eg[hlevelindex];
      normy_ee = G4UniformRand()*_maxWTheta_ee[hlevelindex];

/*	G4cout << "_maxWTheta_gg = " << _maxWTheta_gg[hlevelindex] << " normy_gg = " << normy_gg << G4endl;
	G4cout << "_maxWTheta_ge = " << _maxWTheta_ge[hlevelindex] << " normy_ge = " << normy_ge << G4endl;
	G4cout << "_maxWTheta_eg = " << _maxWTheta_eg[hlevelindex] << " normy_eg = " << normy_eg << G4endl;
	G4cout << "_maxWTheta_ee = " << _maxWTheta_ee[hlevelindex] << " normy_ee = " << normy_ee << G4endl;

	G4cout << "---------- gg vector --------------------------------------" << G4endl;
	G4cout << "a2 to a10: " << _a2_gg[hlevelindex] << " " <<  _a4_gg[hlevelindex] << " " << _a6_gg[hlevelindex] << " " << _a8_gg[hlevelindex] << " " << _a10_gg[hlevelindex] << G4endl;

	G4cout << "---------- ge vector --------------------------------------" << G4endl;
	G4cout << "a2 to a10: " << _a2_ge[hlevelindex] << " " << _a4_ge[hlevelindex] << " " << _a6_ge[hlevelindex] << " " <<  _a8_ge[hlevelindex] << " " << _a10_ge[hlevelindex] << G4endl;

	G4cout << "---------- eg vector --------------------------------------" << G4endl;
	G4cout << "a2 to a10: " << _a2_eg[hlevelindex] << " " <<  _a4_eg[hlevelindex] << " " << _a6_eg[hlevelindex] << " " << _a8_eg[hlevelindex] << " " << _a10_eg[hlevelindex] << G4endl;

	G4cout << "---------- ee vector --------------------------------------" << G4endl;
	G4cout << "a2 to a10: " << _a2_ee[hlevelindex] << " " << _a4_ee[hlevelindex] << " " <<  _a6_ee[hlevelindex] << " " << _a8_ee[hlevelindex] << " " <<  _a10_ee[hlevelindex] << G4endl;
*/

	// casenumber legend
	// 0 == gamma-gamma
	// 1 == gamma-electron
	// 2 == electron-gamma
	// 3 == electron-electron

	// Rishita TODO: add three other cases -- are other thetas necessary?
	
      if (casenumber == 0) {
//	G4cout << "G4DiscreteGammaTransition: case 0" << G4endl;
      	WTheta_gg = (1 + _a2_gg[hlevelindex] * LegendreP(2,std::cos(theta)) 
		  + _a4_gg[hlevelindex] * LegendreP(4,std::cos(theta)) 
 	 	  + _a6_gg[hlevelindex] * LegendreP(6,std::cos(theta))
		  + _a8_gg[hlevelindex] * LegendreP(8,std::cos(theta))
		  + _a10_gg[hlevelindex] * LegendreP(10,std::cos(theta))) * std::sin(theta);

		if(_a2_gg[hlevelindex] == 0 || isnan(_a2_gg[hlevelindex])) { break; } // Rishita -- break or continue? value?

//		value = (WTheta_gg) / _maxWTheta_gg[hlevelindex];
//		break;
	
		if(normy_gg <= WTheta_gg) {
        		value = theta;
        		foundTheta = true;
	//     		G4cout << "CASE 0: WTheta = " << WTheta_gg << " AND normy = " << normy_gg << G4endl;
	//		G4cout << "a2 used: " << _a2_gg[hlevelindex] << G4endl;	
        		break;
      		}
   
      } else if (casenumber == 1) {	
//	G4cout << "G4DiscreteGammaTransition: case 1" << G4endl;
      	WTheta_ge = (1 + _a2_ge[hlevelindex] * LegendreP(2,std::cos(theta)) 
		  + _a4_ge[hlevelindex] * LegendreP(4,std::cos(theta)) 
 	 	  + _a6_ge[hlevelindex] * LegendreP(6,std::cos(theta))
		  + _a8_ge[hlevelindex] * LegendreP(8,std::cos(theta))
		  + _a10_ge[hlevelindex] * LegendreP(10,std::cos(theta))) * std::sin(theta);

		if(_a2_ge[hlevelindex] == 0 || isnan(_a2_ge[hlevelindex])) { break; } // Rishita -- break or continue? value?
		

//		value = (WTheta_ge) / _maxWTheta_ge[hlevelindex];
//		break;

		if(normy_ge <= WTheta_ge) {
        		value = theta;
        		foundTheta = true;
  //     			G4cout << "CASE 1: WTheta = " << WTheta_ge << " AND normy = " << normy_ge << G4endl;	
        		break;
      		}
    
      } else if (casenumber == 2) {
//	G4cout << "G4DiscreteGammaTransition: case 2" << G4endl;
      	WTheta_eg = (1 + _a2_eg[hlevelindex] * LegendreP(2,std::cos(theta)) 
		  + _a4_eg[hlevelindex] * LegendreP(4,std::cos(theta)) 
 	 	  + _a6_eg[hlevelindex] * LegendreP(6,std::cos(theta))
		  + _a8_eg[hlevelindex] * LegendreP(8,std::cos(theta))
		  + _a10_eg[hlevelindex] * LegendreP(10,std::cos(theta))) * std::sin(theta);
      
		if(_a2_eg[hlevelindex] == 0 || isnan(_a2_eg[hlevelindex])) { break; } // Rishita -- break or continue? value?	

//		value = (WTheta_eg) / _maxWTheta_eg[hlevelindex];
//		break;

		if(normy_eg <= WTheta_eg) {
        		value = theta;
        		foundTheta = true;
  //     			G4cout << "CASE 2: WTheta = " << WTheta_eg << " AND normy = " << normy_eg << G4endl;	
        		break;
      		}

      } else if (casenumber == 3) {
//	G4cout << "G4DiscreteGammaTransition: case 3" << G4endl;
      	WTheta_ee = (1 + _a2_ee[hlevelindex] * LegendreP(2,std::cos(theta)) 
		  + _a4_ee[hlevelindex] * LegendreP(4,std::cos(theta)) 
 	 	  + _a6_ee[hlevelindex] * LegendreP(6,std::cos(theta))
		  + _a8_ee[hlevelindex] * LegendreP(8,std::cos(theta))
     		  + _a10_ee[hlevelindex] * LegendreP(10,std::cos(theta))) * std::sin(theta);
      
		if(_a2_ee[hlevelindex] == 0 || isnan(_a2_ee[hlevelindex])) { break; } // Rishita -- break or continue? value?
	
//		value = (WTheta_ee) / _maxWTheta_ee[hlevelindex];
//		break;

		if(normy_ee <= WTheta_ee) {
        		value = theta;
        		foundTheta = true;
//			G4cout << "CASE 3: WTheta = " << WTheta_ee << " AND normy = " << normy_ee << G4endl;	
        		break;
      		}
      }

    }
  }

//G4cout << "------------ G4DiscreteGammaTransition END -------------------------------------------------" << G4endl;


  return value;
}

G4double G4DiscreteGammaTransition::LegendreP(G4int n, G4double x)
{
  G4double PP;
  if(n == 0) {
    PP = 1;
  }
  else if(n == 2) {
    PP = (1.0/2.0)*(3*(std::pow(x,2))-1);
  }
  else if(n == 4) {
    PP = (1.0/8.0)*(35*(std::pow(x,4))-30*(std::pow(x,2))+3);
  }
  else if(n == 6) {
    PP = (1.0/16.0)*(231*(std::pow(x,6))-315*(std::pow(x,4))+105*(std::pow(x,2))-5);
  }
  else if(n == 8) {
    PP = (1.0/128.0)*(6435*(std::pow(x,8))-12012*(std::pow(x,6))+6930*(std::pow(x,4))-1260*(std::pow(x,2))+35);
  }
  else if(n == 10) {
    PP = (1.0/256.0)*(46189*(std::pow(x,10))-109395*(std::pow(x,8))+90090*(std::pow(x,6))-30030*(std::pow(x,4))+3465*(std::pow(x,2))-63);
  }
  else {
    G4cout << "Legendre Polynomial Not Found" << G4endl;
    exit(1);
  }
  return PP;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
