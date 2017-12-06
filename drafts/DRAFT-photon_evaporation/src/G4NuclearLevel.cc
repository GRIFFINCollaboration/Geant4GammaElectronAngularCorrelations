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
// $Id: G4NuclearLevel.cc 86986 2014-11-21 13:00:05Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4NuclearLevel
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//
//      Creation date: 24 October 1998
//
//      Modifications:
//	  06 Oct 2010, M. Kelsey (kelsey@slac.stanford.edu)
//		Add friendship for G4NuclearLevelManager; define private
//		constructors without vectors.
//
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added IC probability when calculate the channel probabilities in
//              MakeProbabilities().
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions.
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//
//        28 October 2010, V.Ivanchenko moved copy constructor to source, cleanup
//
// -------------------------------------------------------------------

#include "G4NuclearLevel.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclearLevelStore.hh" // Will Ashfield 




G4int G4NuclearLevel::Increment(G4int aF)
{
  static G4ThreadLocal G4int instanceCount = 0;
  instanceCount+=aF;
  return instanceCount;
}

G4NuclearLevel::G4NuclearLevel()
  : _energy(0.), _halfLife(0.), _angularMomentum(0.), _nGammas(0) {
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::G4NuclearLevel(G4double energy, G4double halfLife,
                   G4double angularMomentum)
  : _energy(energy), _halfLife(halfLife), _angularMomentum(angularMomentum),
    _nGammas(0) {
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::G4NuclearLevel(G4double energy, G4double halfLife,
                   G4double angularMomentum,
                   const std::vector<G4double>& eGamma,
                   const std::vector<G4double>& wGamma,
                   const std::vector<G4double>& polarities,
                   const std::vector<G4double>& kCC, const std::vector<G4double>& l1CC,
                   const std::vector<G4double>& l2CC, const std::vector<G4double>& l3CC,
                   const std::vector<G4double>& m1CC, const std::vector<G4double>& m2CC,
                   const std::vector<G4double>& m3CC, const std::vector<G4double>& m4CC,
                   const std::vector<G4double>& m5CC, const std::vector<G4double>& nPlusCC,
                   const std::vector<G4double>& totalCC)

  : _energies(eGamma), _weights(wGamma), _polarities(polarities),
     _kCC(kCC), _l1CC(l1CC), _l2CC(l2CC), _l3CC(l3CC),
    _m1CC(m1CC), _m2CC(m2CC), _m3CC(m3CC), _m4CC(m4CC), _m5CC(m5CC),
    _nPlusCC(nPlusCC), _totalCC(totalCC),
    _energy(energy), _halfLife(halfLife), _angularMomentum(angularMomentum)
{
  Finalize();
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::~G4NuclearLevel()
{
 // G4cout << "####### Decrementing "<<Increment(-1)<<G4endl;
}

G4bool G4NuclearLevel::operator==(const G4NuclearLevel &right) const
{
  return (this == (G4NuclearLevel *) &right);
}


G4bool G4NuclearLevel::operator!=(const G4NuclearLevel &right) const
{
  return (this != (G4NuclearLevel *) &right);
}

G4bool G4NuclearLevel::operator<(const G4NuclearLevel &right) const
{
  if (_energy < right.Energy()) return true;
  else return false;
}

const std::vector<G4double>& G4NuclearLevel::GammaEnergies() const
{
  return _energies;
}

const std::vector<G4double>& G4NuclearLevel::GammaWeights() const
{
  return _weights;
}


const std::vector<G4double>& G4NuclearLevel::GammaProbabilities() const
{
  return _prob;
}


const std::vector<G4double>& G4NuclearLevel::GammaCumulativeProbabilities() const
{
  return _cumProb;
}


const std::vector<G4double>& G4NuclearLevel::GammaPolarities() const
{
  return _polarities;
}

const std::vector<G4double>& G4NuclearLevel::KConvertionProbabilities() const
{
  return _kCC;
}

const std::vector<G4double>& G4NuclearLevel::L1ConvertionProbabilities() const
{
  return _l1CC;
}

const std::vector<G4double>& G4NuclearLevel::L2ConvertionProbabilities() const
{
  return _l2CC;
}

const std::vector<G4double>& G4NuclearLevel::L3ConvertionProbabilities() const
{
  return _l3CC;
}

const std::vector<G4double>& G4NuclearLevel::M1ConvertionProbabilities() const
{
  return _m1CC;
}

const std::vector<G4double>& G4NuclearLevel::M2ConvertionProbabilities() const
{
  return _m2CC;
}

const std::vector<G4double>& G4NuclearLevel::M3ConvertionProbabilities() const
{
  return _m3CC;
}

const std::vector<G4double>& G4NuclearLevel::M4ConvertionProbabilities() const
{
  return _m4CC;
}

const std::vector<G4double>& G4NuclearLevel::M5ConvertionProbabilities() const
{
  return _m5CC;
}

const std::vector<G4double>& G4NuclearLevel::NPlusConvertionProbabilities() const
{
  return _nPlusCC;
}

const std::vector<G4double>& G4NuclearLevel::TotalConvertionProbabilities() const
{
  return _totalCC;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

const std::vector<G4double>& G4NuclearLevel::L1() const
{
  return _L1;
}

const std::vector<G4double>& G4NuclearLevel::L2() const
{
  return _L2;
}

const std::vector<G4double>& G4NuclearLevel::MixingRatio() const
{
  return _mixingRatio;
}

// Rishita ////////////////////////////////////////////////////////////////////////////////////////////

const std::vector<G4double>& G4NuclearLevel::B1() const
{
  return _b1;
}

const std::vector<G4double>& G4NuclearLevel::B2() const
{
  return _b2;
}

const std::vector<G4double>& G4NuclearLevel::B3() const
{
  return _b3;
}

const std::vector<G4double>& G4NuclearLevel::B4() const
{
  return _b4;
}

const std::vector<G4double>& G4NuclearLevel::B5() const
{
  return _b5;
}

const std::vector<G4double>& G4NuclearLevel::B6() const
{
  return _b6;
}

const std::vector<G4double>& G4NuclearLevel::AlphL1() const
{
  return _alphL1;
}

const std::vector<G4double>& G4NuclearLevel::AlphL2() const
{
  return _alphL2;
}

const std::vector<G4double>& G4NuclearLevel::Electrue() const
{
  return _electrue;
} 

// gamma-gamma -----------------------------------------------------------------------------------------
const std::vector< std::vector< G4double > >& G4NuclearLevel::A2_GG() const
{
  return _a2_gg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A4_GG() const
{
  return _a4_gg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A6_GG() const
{
  return _a6_gg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A8_GG() const
{
  return _a8_gg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A10_GG() const
{
  return _a10_gg;
}

// gamma-electron --------------------------------------------------------------------------------------
const std::vector< std::vector< G4double > >& G4NuclearLevel::A2_GE() const
{
  return _a2_ge;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A4_GE() const
{
  return _a4_ge;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A6_GE() const
{
  return _a6_ge;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A8_GE() const
{
  return _a8_ge;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A10_GE() const
{
  return _a10_ge;
}

// electron-gamma --------------------------------------------------------------------------------------
const std::vector< std::vector< G4double > >& G4NuclearLevel::A2_EG() const
{
  return _a2_eg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A4_EG() const
{
  return _a4_eg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A6_EG() const
{
  return _a6_eg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A8_EG() const
{
  return _a8_eg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A10_EG() const
{
  return _a10_eg;
}

// electron-electron -----------------------------------------------------------------------------------
const std::vector< std::vector< G4double > >& G4NuclearLevel::A2_EE() const
{
  return _a2_ee;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A4_EE() const
{
  return _a4_ee;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A6_EE() const
{
  return _a6_ee;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A8_EE() const
{
  return _a8_ee;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::A10_EE() const
{
  return _a10_ee;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////

const std::vector< std::vector< G4double > >& G4NuclearLevel::MaxWTheta_GG() const
{
  return _maxWTheta_gg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::MaxWTheta_GE() const
{
  return _maxWTheta_ge;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::MaxWTheta_EG() const
{
  return _maxWTheta_eg;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::MaxWTheta_EE() const
{
  return _maxWTheta_ee;
}

const std::vector< std::vector< G4double > >& G4NuclearLevel::HigherLevelEnergy() const
{
  return _higherLevelEnergy;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

G4double G4NuclearLevel::Energy() const
{
  return _energy;
}

G4double G4NuclearLevel::AngularMomentum() const
{
  return _angularMomentum;
}

G4double G4NuclearLevel::HalfLife() const
{
  return _halfLife;
}

G4int G4NuclearLevel::NumberOfGammas() const
{
  return _nGammas;
}

void G4NuclearLevel::PrintAll() const
{
  G4cout << "---- Level energy = " << _energy << ", angular momentum = "
     << _angularMomentum << ", half life " << _halfLife
     << ", " << _nGammas << " photons" << G4endl;
  G4int i;
  G4cout << "     Gammas: ";
  for (i=0; i<_nGammas; i++) { G4cout << _energies[i] << " "; }
  G4cout << G4endl << "     Weights: ";
  for (i=0; i<_nGammas; i++) { G4cout << _weights[i] << " "; }
  G4cout << G4endl << "     Relative transition probabilities ";
  for (i=0; i<_nGammas; i++) { G4cout << _prob[i] << " "; }
  G4cout << G4endl << "     Cumulative probabilities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _cumProb[i] << " "; }
  G4cout << G4endl << "     Polarities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _polarities[i] << " "; }
  G4cout << G4endl;
}

void G4NuclearLevel::PrintLevels() const
{
  G4cout << "   Eexc(MeV)= " << _energy
     << " Time(ns)= " << _halfLife/ns << "  Ntrans= " << _nGammas
     << G4endl;
}

void G4NuclearLevel::Finalize() {
  _nGammas = _energies.size();
  MakeProbabilities();
  MakeCumProb();
  MakeMultipoleData(); // Evan Rand
}

void G4NuclearLevel::MakeProbabilities()
{
  G4double sum = 0.;
  G4int i = 0;
  for (i=0; i<_nGammas; i++) {
    sum += _weights[i]*(1.+_totalCC[i]);
  }

  if (sum <= 0.) _prob.resize(_nGammas, 1./_nGammas);	// Fast fill
  else {
    _prob.reserve(_nGammas);
    for (i=0; i<_nGammas; i++) {
      _prob.push_back(_weights[i]*(1.+_totalCC[i])/sum);
    }
  }
}


void G4NuclearLevel::MakeCumProb()
{
  if (_nGammas <= 0) return;

  _cumProb.reserve(_nGammas);

  G4double sum = _prob[0];
  _cumProb.push_back(sum);

  for (G4int i=1; i<_nGammas; i++) {
    sum += _prob[i];
    _cumProb.push_back(sum);
  }
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void G4NuclearLevel::MakeMultipoleData()
{
  if (_nGammas <= 0) return;

  for (G4int i=0; i<_nGammas; i++) {
    _L1.push_back(0);
    _L2.push_back(0);
    _mixingRatio.push_back(0);

// Rishita --------------------------------------------------------------------------------------------
	
    _b1.push_back(0);
    _b2.push_back(0);
    _b3.push_back(0);
    _b4.push_back(0);
    _b5.push_back(0);
    _b6.push_back(0);
    _alphL1.push_back(0);
    _alphL2.push_back(0);
    _electrue.push_back(0); 
   
	// TODO: add vectors for other three cases 
    _higherLevelEnergy.push_back(std::vector< G4double >());
   
    _a2_gg.push_back(std::vector< G4double >());
    _a4_gg.push_back(std::vector< G4double >());
    _a6_gg.push_back(std::vector< G4double >());
    _a8_gg.push_back(std::vector< G4double >());
    _a10_gg.push_back(std::vector< G4double >());

    _a2_ge.push_back(std::vector< G4double >());
    _a4_ge.push_back(std::vector< G4double >());
    _a6_ge.push_back(std::vector< G4double >());
    _a8_ge.push_back(std::vector< G4double >());
    _a10_ge.push_back(std::vector< G4double >());
   
    _a2_eg.push_back(std::vector< G4double >());
    _a4_eg.push_back(std::vector< G4double >());
    _a6_eg.push_back(std::vector< G4double >());
    _a8_eg.push_back(std::vector< G4double >());
    _a10_eg.push_back(std::vector< G4double >());
   
    _a2_ee.push_back(std::vector< G4double >());
    _a4_ee.push_back(std::vector< G4double >());
    _a6_ee.push_back(std::vector< G4double >());
    _a8_ee.push_back(std::vector< G4double >());
    _a10_ee.push_back(std::vector< G4double >());
   
    _maxWTheta_gg.push_back(std::vector< G4double >());
    _maxWTheta_ge.push_back(std::vector< G4double >());
    _maxWTheta_eg.push_back(std::vector< G4double >());
    _maxWTheta_ee.push_back(std::vector< G4double >());
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

G4NuclearLevel& G4NuclearLevel::operator=(const G4NuclearLevel &right)
{
  if(this != &right)
    {
      _energies = right._energies;
      _weights =right._weights;
      _prob =right._prob;
      _cumProb =right._cumProb;
      _polarities =right._polarities;
      _kCC = right._kCC;
      _l1CC =right._l1CC;
      _l2CC =right._l2CC;
      _l3CC =right._l3CC;
      _m1CC = right._m1CC;
      _m2CC = right._m2CC;
      _m3CC = right._m3CC;
      _m4CC = right._m4CC;
      _m5CC = right._m5CC;
      _nPlusCC = right._nPlusCC;
      _totalCC = right._totalCC;
      _energy = right._energy;
      _halfLife = right._halfLife;
      _angularMomentum = right._angularMomentum;
      _nGammas = right._nGammas;
    }
  return *this;
}

G4NuclearLevel::G4NuclearLevel(const G4NuclearLevel &right)
{
  _energies = right._energies;
  _weights =right._weights;
  _prob =right._prob;
  _cumProb =right._cumProb;
  _polarities =right._polarities;
  _kCC = right._kCC;
  _l1CC =right._l1CC;
  _l2CC =right._l2CC;
  _l3CC =right._l3CC;
  _m1CC = right._m1CC;
  _m2CC = right._m2CC;
  _m3CC = right._m3CC;
  _m4CC = right._m4CC;
  _m5CC = right._m5CC;
  _nPlusCC = right._nPlusCC;
  _totalCC = right._totalCC;
  _energy = right._energy;
  _halfLife = right._halfLife;
  _angularMomentum = right._angularMomentum;
  _nGammas = right._nGammas;
}

// Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void G4NuclearLevel::FillHigherLevelEnergy(G4int gamma_i, G4int hlevel_i, G4double hlevel_energy)
{
  G4int nLevels = _higherLevelEnergy[gamma_i].size();
  if (hlevel_i > nLevels){
    G4cout << "Error! Skipped a level?" << G4endl;
    exit(1);
  }
  _higherLevelEnergy[gamma_i].push_back(hlevel_energy);
}

void G4NuclearLevel::GenerateWThetaParameters(G4int gamma_i, G4int hlevel_i, G4double hlevel_energy, G4double lowerGammaEnergy, G4double higherGammaEnergy, G4double ji, G4double jo, G4double jf, G4int L1, G4int L1p, G4int L2, G4int L2p, G4double delta1, G4double delta2, G4bool boolGoodLevelToOutputToScreen, G4double b_a1, G4double b_a2, G4double b_a3, G4double b_a4, G4double b_a5, G4double b_a6, G4double alph_aL, G4double alph_aLp, G4double electrue_a, G4double b_b1, G4double b_b2, G4double b_b3, G4double b_b4, G4double b_b5, G4double b_b6, G4double alph_bL, G4double alph_bLp, G4double electrue_b) {
 
  G4double P2, P4, P6, P8, P10, thetaRad;
  G4double a2_gg, a4_gg, a6_gg, a8_gg, a10_gg, wTheta_gg, thisMaxWTheta_gg;
  G4double a2_ge, a4_ge, a6_ge, a8_ge, a10_ge, wTheta_ge, thisMaxWTheta_ge;
  G4double a2_eg, a4_eg, a6_eg, a8_eg, a10_eg, wTheta_eg, thisMaxWTheta_eg;
  G4double a2_ee, a4_ee, a6_ee, a8_ee, a10_ee, wTheta_ee, thisMaxWTheta_ee;

  // a check for non-sensical spin/mixing ratio combinations -SmithJK
  if (L1==L1p) delta1 = 0;
  if (L2==L2p) delta2 = 0;

  // Save higher level energy
  FillHigherLevelEnergy(gamma_i, hlevel_i, hlevel_energy);

  // An additional check that these variables are positive
  G4double Ji = std::fabs(ji);
  G4double Jo = std::fabs(jo);
  G4double Jf = std::fabs(jf);

  // Rishita ---------------------------------------------------------------------

  G4double delta_ex1 = std::pow((alph_aLp / alph_aL), 0.5)*delta1;
  G4double delta_ex2 = std::pow((alph_bLp / alph_bL), 0.5)*delta2;

  // gamma-gamma transition
  a2_gg  = firstTransition(2,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(2,Jf,Jo,L2,L2p,delta2,1,1,1);
  a4_gg  = firstTransition(4,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(4,Jf,Jo,L2,L2p,delta2,1,1,1);
  a6_gg  = firstTransition(6,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(6,Jf,Jo,L2,L2p,delta2,1,1,1);
  a8_gg  = firstTransition(8,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(8,Jf,Jo,L2,L2p,delta2,1,1,1);
  a10_gg = firstTransition(10,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(10,Jf,Jo,L2,L2p,delta2,1,1,1);

	G4double temp1 = firstTransition(2,Jo,Ji,L1,L1p,delta1,1,1,1);
	G4double temp2 = secondTransition(2,Jf,Jo,L2,L2p,delta2,1,1,1);
//	G4cout << "---------------------------------------------------------" << G4endl;
//	G4cout << "a2_gg firstTransition: " << temp1 << G4endl;
//	G4cout << "Jo: " << Jo << ", Ji: " << Ji << ", L1: " << L1 << ", L1p: " << L1p << ", delta1: " << delta1 << G4endl;
//	G4cout << "a2_gg secondTransition: " << temp2 << G4endl;
//	G4cout << "Jf: " << Jf << ", Jo: " << Jo << ", L2: " << L2 << ", L2p: " << L2p << ", delta2: " << delta2 << G4endl;
//	G4cout << "------------------------------------------------------" << G4endl;

  // gamma-electron transition
  a2_ge  = firstTransition(2,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(2,Jf,Jo,L2,L2p,delta_ex2,b_b1,b_b2,b_b3);
  a4_ge  = firstTransition(4,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(4,Jf,Jo,L2,L2p,delta_ex2,b_b4,b_b5,b_b6);
  a6_ge  = firstTransition(6,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(6,Jf,Jo,L2,L2p,delta_ex2,1,1,1);
  a8_ge  = firstTransition(8,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(8,Jf,Jo,L2,L2p,delta_ex2,1,1,1);
  a10_ge = firstTransition(10,Jo,Ji,L1,L1p,delta1,1,1,1)*secondTransition(10,Jf,Jo,L2,L2p,delta_ex2,1,1,1);

  // electron-gamma transition
  a2_eg  = firstTransition(2,Jo,Ji,L1,L1p,delta_ex1,b_a1,b_a2,b_a3)*secondTransition(2,Jf,Jo,L2,L2p,delta2,1,1,1);
  a4_eg  = firstTransition(4,Jo,Ji,L1,L1p,delta_ex1,b_a4,b_a5,b_a6)*secondTransition(4,Jf,Jo,L2,L2p,delta2,1,1,1);
  a6_eg  = firstTransition(6,Jo,Ji,L1,L1p,delta_ex1,1,1,1)*secondTransition(6,Jf,Jo,L2,L2p,delta2,1,1,1);
  a8_eg  = firstTransition(8,Jo,Ji,L1,L1p,delta_ex1,1,1,1)*secondTransition(8,Jf,Jo,L2,L2p,delta2,1,1,1);
  a10_eg = firstTransition(10,Jo,Ji,L1,L1p,delta_ex1,1,1,1)*secondTransition(10,Jf,Jo,L2,L2p,delta2,1,1,1);
  
  // electron-electron transition
  a2_ee  = firstTransition(2,Jo,Ji,L1,L1p,delta_ex1,b_a1,b_a2,b_a3)*secondTransition(2,Jf,Jo,L2,L2p,delta_ex2,b_b1,b_b2,b_b3);
  a4_ee  = firstTransition(4,Jo,Ji,L1,L1p,delta_ex1,b_a4,b_a5,b_a6)*secondTransition(4,Jf,Jo,L2,L2p,delta_ex2,b_b4,b_b5,b_b6);
  a6_ee  = firstTransition(6,Jo,Ji,L1,L1p,delta_ex1,1,1,1)*secondTransition(6,Jf,Jo,L2,L2p,delta_ex2,1,1,1);
  a8_ee  = firstTransition(8,Jo,Ji,L1,L1p,delta_ex1,1,1,1)*secondTransition(8,Jf,Jo,L2,L2p,delta_ex2,1,1,1);
  a10_ee = firstTransition(10,Jo,Ji,L1,L1p,delta_ex1,1,1,1)*secondTransition(10,Jf,Jo,L2,L2p,delta_ex2,1,1,1);
  

  // Will Ashfield - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Rishita TODO - set manual coefficients for other three cases
  if(G4NuclearLevelStore::GetInstance()->manualACcoeffs()){
    a2_gg = G4NuclearLevelStore::GetInstance()->GetA2();  
    a4_gg = G4NuclearLevelStore::GetInstance()->GetA4(); 
    a6_gg = G4NuclearLevelStore::GetInstance()->GetA6();
    a8_gg = 0;
    a10_gg = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  

  if(boolGoodLevelToOutputToScreen && (a2_gg != 0 || a4_gg != 0)) {
    G4cout << "------------------- angular coefficients --------------------" << G4endl;
    G4cout << "highest level in cascade is " << hlevel_energy/keV << " keV" << G4endl;
    G4cout << "first gamma in cascade is " << higherGammaEnergy/keV << " keV" << G4endl;
    G4cout << "second gamma in cascade is " << lowerGammaEnergy/keV << " keV" << G4endl;
    if(G4NuclearLevelStore::GetInstance()->manualACcoeffs()){
    G4cout << "THE ANGULAR COEFFICIENTS HAVE BEEN MANUALLY SET" << G4endl;  
    G4cout << "The simulation will now overwrite the calculated a2 and a4 coefficients" << G4endl;  
    G4cout << "This process will overwrite the coefficients for any cascade" << G4endl;  
    } else{
    G4cout << "ji = " << ji << " --> jo = " << jo <<  " --> jf = " << jf << G4endl;
    G4cout << "L1 = " << L1 << ", L1p = " << L1p <<  " - " << "L2 = " << L2 << ", L2p = " << L2p << G4endl;
    G4cout << "delta1 (gamma) = " << delta1 << " - " << "delta2 (gamma) = " << delta2 << G4endl;
    G4cout << "delta1 (electron) = " << delta_ex1 << " - " << "delta2 (electron) = " << delta_ex2 << G4endl;
    }
    G4cout << "------------------- gamma-gamma transition ------------------" << G4endl;
    G4cout << "a2  = " << a2_gg << G4endl;
    G4cout << "a4  = " << a4_gg << G4endl;
    G4cout << "a6  = " << a6_gg << G4endl;
    G4cout << "a8  = " << a8_gg << G4endl;
    G4cout << "a10 = " << a10_gg << G4endl;
    G4cout << "------------------- gamma-electron transition ---------------" << G4endl;
    G4cout << "a2  = " << a2_ge << G4endl;
    G4cout << "a4  = " << a4_ge << G4endl;
    G4cout << "a6  = " << a6_ge << G4endl;
    G4cout << "a8  = " << a8_ge << G4endl;
    G4cout << "a10 = " << a10_ge << G4endl;
    G4cout << "------------------ electron-gamma transition ----------------" << G4endl;
    G4cout << "a2  = " << a2_eg << G4endl;
    G4cout << "a4  = " << a4_eg << G4endl;
    G4cout << "a6  = " << a6_eg << G4endl;
    G4cout << "a8  = " << a8_eg << G4endl;
    G4cout << "a10 = " << a10_eg << G4endl;
    G4cout << "------------------ electron-electron transition -------------" << G4endl;
    G4cout << "a2  = " << a2_ee << G4endl;
    G4cout << "a4  = " << a4_ee << G4endl;
    G4cout << "a6  = " << a6_ee << G4endl;
    G4cout << "a8  = " << a8_ee << G4endl;
    G4cout << "a10 = " << a10_ee << G4endl;
    G4cout << "-------------------------------------------------------------" << G4endl;
    if(G4NuclearLevelStore::GetInstance()->manualACcoeffs()){
    G4cout << "REMINDER: THE ANGULAR COEFFICIENTS HAVE BEEN MANUALLY SET" << G4endl;
    }/*
    G4cout << "-------------------------------- gamma-electron angular coefficients -------------------------------------" << G4endl;
    G4cout << "Higher Energy Parameters are found on 2nd line of Multipole File & have subscript a for first transition" << G4endl;
    G4cout << "Lower Energy Parameters are found on 1st line of Multipole File & have subscript b for second transition" << G4endl;
    G4cout << "----------------------------------------------------------------------------------------------------------" << G4endl;
    G4cout << ".........FIRST TRANSITION........." << G4endl;
    G4cout << "b1 to b6 passed: " << b_a1 << ", " << b_a2 << ", " << b_a3 << ", " << b_a4 << ", " << b_a5 << ", " << b_a6 << "." << G4endl;
    G4cout << "b1 to b6 used: " << b1a_val << ", " << b2a_val << ", " << b3a_val << ", " << b4a_val << ", " << b5a_val << ", " << b6a_val << "." <<G4endl;
    G4cout << "alpha_L: " <<  alph_aL << G4endl;
    G4cout << "alpha_Lp: " << alph_aLp << G4endl;
    	if (ea) { G4cout << "*****first transition was an e-*****" << G4endl; }
    	else { G4cout << "*****first transition was a gamma*****" << G4endl; }
    G4cout << "------------------------------------------------------" << G4endl;
    G4cout << ".........SECOND TRANSITION........." << G4endl;
    G4cout << "b1 to b6 passed: " << b_b1 << ", " << b_b2 << ", " << b_b3 << ", " << b_b4 << ", " << b_b5 << ", " << b_b6 << "." << G4endl;
    G4cout << "trans b1 to b6 used: " << b1b_val << ", " << b2b_val << ", " << b3b_val << ", " << b4b_val << ", " << b5b_val << ", " << b6b_val << "." <<G4endl;
    G4cout << "alpha_L: " <<  alph_bL << G4endl;
    G4cout << "alpha_Lp: " << alph_bLp << G4endl;
    	if (eb) { G4cout << "*****second transition was an e-*****" << G4endl; }
    	else { G4cout << "*****second transition was a gamma*****" << G4endl; }
    G4cout << "------------------------------------------------------" << G4endl;
    G4cout << "electrue_a: " << electrue_a << G4endl;
    G4cout << "electrue_b: " << electrue_b << G4endl;
    G4cout << "------------------------------------------------------" << G4endl;*/
    //std::vector<const G4ParticleDefinition*>* rep = G4NuclearLevelStore::GetInstance()->GetResults();
    /*if(!rep) {G4cout <<"nothing"<<G4endl;} 
	else {
    for(G4int i = 0; i < rep->size(); i++) {
        const G4ParticleDefinition* n = rep->at(i);
        const G4String name = n->GetParticleName();
	G4cout << "Particle Name: " << name << G4endl;
    }}*/
    // DELETE THIS ----------------------------------------------------------------
    /*G4ParticleTable::G4PTblDicIterator* theParticleIterator
        = G4ParticleTable::GetParticleTable()->GetIterator();

    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
        G4ParticleDefinition* particle = theParticleIterator->value();
	G4String particleName = particle->GetParticleName();
        G4cout << "particleName: " << particleName << G4endl;
    }*/
    //G4ParticleDefinition* particle = G4ParticleGun::GetParticleDefinition();
    //G4cout << "particleName: " << particle->GetParticleName() << G4endl;
//-----------------------------------------------------------------------------

//const G4String Def0 = ((*secondaries)[0])->GetDefinition()->GetParticleName();
//const G4String Def1 = ((*secondaries)[1])->GetDefinition()->GetParticleName();
   // G4cout << "Def0: " << Def0 << G4endl;
    //G4cout << "Def1: " << Def1 << G4endl;

    /*if (caseno == 0) { G4cout << "no case assigned" << G4endl; }
    else if (caseno == 1) { G4cout << "electrue_a & electrue_b both false" << G4endl; }
    else if (caseno == 2) { G4cout << "electrue_a false & electrue_b true" << G4endl; }
    else if (caseno == 3) { G4cout << "electrue_a true & electrue_b false" << G4endl; }
    else if (caseno == 4) { G4cout << "electrue_a & electrue_b both true" << G4endl; }
    else { G4cout << "???" << G4endl; }
    G4cout << "-------------------------------------------------------------------------------------------" << G4endl;*/
}

  // Rishita -- TODO: add vectors for other three cases
  _a2_gg[gamma_i].push_back(a2_gg);
  _a4_gg[gamma_i].push_back(a4_gg);
  _a6_gg[gamma_i].push_back(a6_gg);
  _a8_gg[gamma_i].push_back(a8_gg);
  _a10_gg[gamma_i].push_back(a10_gg);
  
  _a2_ge[gamma_i].push_back(a2_ge);
  _a4_ge[gamma_i].push_back(a4_ge);
  _a6_ge[gamma_i].push_back(a6_ge);
  _a8_ge[gamma_i].push_back(a8_ge);
  _a10_ge[gamma_i].push_back(a10_ge);
  
  _a2_eg[gamma_i].push_back(a2_eg);
  _a4_eg[gamma_i].push_back(a4_eg);
  _a6_eg[gamma_i].push_back(a6_eg);
  _a8_eg[gamma_i].push_back(a8_eg);
  _a10_eg[gamma_i].push_back(a10_eg);
  
  _a2_ee[gamma_i].push_back(a2_ee);
  _a4_ee[gamma_i].push_back(a4_ee);
  _a6_ee[gamma_i].push_back(a6_ee);
  _a8_ee[gamma_i].push_back(a8_ee);
  _a10_ee[gamma_i].push_back(a10_ee);


  // find max wTheta to normalize to later.
  // Rishita -- TODO: add max wTheta for other three cases
  thisMaxWTheta_gg = 0;
  thisMaxWTheta_ge = 0;
  thisMaxWTheta_eg = 0;
  thisMaxWTheta_ee = 0;
  for(G4int i = 0; i < 180000 ; i++ ) {
    thetaRad = (i/1000)*(M_PI/180.0);
    P2 = LegendreP(2,cos(thetaRad));
    P4 = LegendreP(4,cos(thetaRad));
    P6 = LegendreP(6,cos(thetaRad));
    P8 = LegendreP(8,cos(thetaRad));
    P10 = LegendreP(10,cos(thetaRad));

    wTheta_gg = (1+a2_gg*P2+a4_gg*P4+a6_gg*P6+a8_gg*P8+a10_gg*P10)*std::sin(thetaRad);
    if (wTheta_gg > thisMaxWTheta_gg) { thisMaxWTheta_gg = wTheta_gg; }

    wTheta_ge = (1+a2_ge*P2+a4_ge*P4+a6_ge*P6+a8_ge*P8+a10_ge*P10)*std::sin(thetaRad);
    if (wTheta_ge > thisMaxWTheta_ge) { thisMaxWTheta_ge = wTheta_ge; }
    
    wTheta_eg = (1+a2_eg*P2+a4_eg*P4+a6_eg*P6+a8_eg*P8+a10_eg*P10)*std::sin(thetaRad);
    if (wTheta_eg > thisMaxWTheta_eg) { thisMaxWTheta_eg = wTheta_eg; }
    
    wTheta_ee = (1+a2_ee*P2+a4_ee*P4+a6_ee*P6+a8_ee*P8+a10_ee*P10)*std::sin(thetaRad);
    if (wTheta_ee > thisMaxWTheta_ee) { thisMaxWTheta_ee = wTheta_ee; }
  }
  // Add a positive 1% error on maxWTheta, this ensures we monte carlo over the entire y-range later.
  // This may be an unnecessary check, but better to be safe than sorry.
  // It will be negligibly slower to add this insurance.

  _maxWTheta_gg[gamma_i].push_back(thisMaxWTheta_gg+(0.01*thisMaxWTheta_gg));
  _maxWTheta_ge[gamma_i].push_back(thisMaxWTheta_ge+(0.01*thisMaxWTheta_ge));
  _maxWTheta_eg[gamma_i].push_back(thisMaxWTheta_eg+(0.01*thisMaxWTheta_eg));
  _maxWTheta_ee[gamma_i].push_back(thisMaxWTheta_ee+(0.01*thisMaxWTheta_ee));
}


G4double G4NuclearLevel::ClebschGordan(G4double j1, G4double m1, G4double j2, G4double m2, G4double j, G4double m)
{
  // Conditions check
  if( 2*j1 !=   std::floor(2*j1) || 2*j2 !=   std::floor(2*j2) || 2*j !=   std::floor(2*j) || 2*m1 !=   std::floor(2*m1) || 2*m2 !=   std::floor(2*m2) || 2*m !=   std::floor(2*m) )
  {
    G4cout << "All arguments must be integers or half-integers." << G4endl;
    return 0;
  }
  if(m1 + m2 != m)
  {
    //G4cout << "m1 + m2 must equal m." << G4endl;
    return 0;
  }
  if( j1 - m1 !=   std::floor ( j1 - m1 ) )
  {
    //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
    return 0;
  }
  if( j2 - m2 !=   std::floor ( j2 - m2 ) )
  {
    //G4cout << "2*j2 and 2*m2 must have the same parity" << G4endl;
    return 0;
  }
  if( j - m !=   std::floor ( j - m ) )
  {
    //G4cout << "2*j and 2*m must have the same parity" << G4endl;
    return 0;
  }
  if(j > j1 + j2 || j < std::fabs(j1 - j2))
  {
    //G4cout << "j is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m1) > j1)
  {
    //G4cout << "m1 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m2) > j2)
  {
    //G4cout << "m2 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m) > j)
  {
    //warning('m is out of bounds." << G4endl;
    return 0 ;
  }
  G4double term, cg;
  G4double term1 = std::pow((((2*j+1)/Factorial(j1+j2+j+1))*Factorial(j2+j-j1)*Factorial(j+j1-j2)*Factorial(j1+j2-j)*Factorial(j1+m1)*Factorial(j1-m1)*Factorial(j2+m2)*Factorial(j2-m2)*Factorial(j+m)*Factorial(j-m)),(0.5));
  G4double sum = 0;

  for(G4int k = 0 ; k <= 99 ; k++ )
  {
    if( (j1+j2-j-k < 0) || (j-j1-m2+k < 0) || (j-j2+m1+k < 0) || (j1-m1-k < 0) || (j2+m2-k < 0) )
    {
    }
    else
    {
      term = Factorial(j1+j2-j-k)*Factorial(j-j1-m2+k)*Factorial(j-j2+m1+k)*Factorial(j1-m1-k)*Factorial(j2+m2-k)*Factorial(k);
      if((k%2) == 1)
      {
        term = -1*term;
      }
      sum = sum + 1.0/term;
    }
  }

  cg = term1*sum;
  return cg;
  // Reference: An Effective Algorithm for Calculation of the C.G.
  // Coefficients Liang Zuo, et. al.
  // J. Appl. Cryst. (1993). 26, 302-304
}

G4double G4NuclearLevel::Factorial(G4double value)
{
  G4double fac;
  if(value > 1)
  {
    fac = value*Factorial(value-1);
  }
  else
  {
    fac = 1;
  }
  return fac;
}


G4double G4NuclearLevel::LegendreP(G4int n, G4double x)
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

G4double G4NuclearLevel::F(G4int k, G4double jf, G4int L1, G4int L2, G4double ji)
{
  G4double out;
  G4double CG = ClebschGordan(L1,1,L2,-1,k,0);
  if(CG == 0)
  {
    return 0;
  }
  G4double W = RacahW(ji,ji,L1,L2,k,jf);
  if(W == 0)
  {
    return 0;
  }
  out = std::pow((-1),(jf-ji-1))*(std::pow((2*L1+1)*(2*L2+1)*(2*ji+1),(1.0/2.0)))*CG*W;
  return out;
  // Reference: Tables of coefficients for angular distribution of gamma rays from aligned nuclei
  // T. Yamazaki. Nuclear Data A, 3(1):1?23, 1967.
}


// Rishita -- modified to include b parameters --> set to 1 for a gamma transition
G4double G4NuclearLevel::firstTransition(G4int k, G4double ji, G4double jf, G4int L1, G4int L2, G4double delta, G4double b1, G4double b2, G4double b3)
{
  G4double out;
  out = (1/(1+std::pow(delta,2)))*(b1* F(k,jf,L1,L1,ji)+(std::pow((-1),((G4double)(L1-L2))))*2*delta*b2*F(k,jf,L1,L2,ji)+delta*delta*b3*F(k,jf,L2,L2,ji) );
  return out;
}


G4double G4NuclearLevel::RacahW(G4double a, G4double b, G4double c, G4double d, G4double e, G4double f)
{
  G4double out = std::pow((-1),(a+b+d+c))*Wigner6j(a,b,e,d,c,f);
  return out;
}

G4double G4NuclearLevel::Wigner6j(G4double J1, G4double J2, G4double J3, G4double J4, G4double J5, G4double J6)
{
  // Conditions check
  if(J3 > J1 + J2 || J3 < std::fabs(J1 - J2))
  {
    //G4cout << "first J3 triange condition not satisfied. J3 > J1 + J2 || J3 < std::fabs(J1 - J2)" << G4endl;
    return 0;
  }
  if(J3 > J4 + J5 || J3 < std::fabs(J4 - J5))
  {
    //G4cout << "second J3 triange condition not satisfied. J3 > J4 + J5 || J3 < std::fabs(J4 - J5)" << G4endl;
    return 0;
  }
  if(J6 > J2 + J4 || J6 < std::fabs(J2 - J4))
  {
    //G4cout << "first J6 triange condition not satisfied. J6 > J2 + J4 || J6 < std::fabs(J2 - J4)" << G4endl;
    return 0;
  }
  if(J6 > J1 + J5 || J6 < std::fabs(J1 - J5))
  {
    //G4cout << "second J6 triange condition not satisfied. J6 > J1 + J5 || J6 < std::fabs(J1 - J5)" << G4endl;
    return 0;
  }

  G4double j1 = J1;
  G4double j2 = J2;
  G4double j12 = J3;
  G4double j3 = J4;
  G4double j = J5;
  G4double j23 = J6;
  G4double sum = 0;

  for(G4double m1 = -j1 ; m1 <= j1 ; m1++ )
  {
    for(G4double m2 = -j2 ; m2 <= j2 ; m2++ )
    {
      for(G4double m3 = -j3 ; m3 <= j3 ; m3++ )
      {
        for(G4double m12 = -j12 ; m12 <= j12 ; m12++ )
        {
          for(G4double m23 = -j23 ; m23 <= j23 ; m23++ )
          {
            for(G4double m = -j ; m <= j ; m++ )
            {
              sum = sum + std::pow((-1),(j3+j+j23-m3-m-m23))*Wigner3j(j1,j2,j12,m1,m2,m12)*Wigner3j(j1,j,j23,m1,-m,m23)*Wigner3j(j3,j2,j23,m3,m2,-m23)*Wigner3j(j3,j,j12,-m3,m,m12);
            }
          }
        }
      }
    }
  }
  return sum;
}

G4double G4NuclearLevel::Wigner3j(G4double j1, G4double j2, G4double j3, G4double m1, G4double m2, G4double m3)
{
  // Conditions check
  if( 2*j1 !=   std::floor(2*j1) || 2*j2 !=   std::floor(2*j2) || 2*j3 !=   std::floor(2*j3) || 2*m1 !=   std::floor(2*m1) || 2*m2 !=   std::floor(2*m2) || 2*m3 !=   std::floor(2*m3) )
  {
    G4cout << "All arguments must be integers or half-integers." << G4endl;
    return 0;
  }
  if(m1 + m2 + m3 != 0)
  {
    //G4cout << "m1 + m2 + m3 must equal zero." << G4endl;
    return 0;
  }
  if( j1 + j2 + j3 !=   std::floor(j1 + j2 + j3) )
  {
    //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
    return 0;
  }
  if(j3 > j1 + j2 || j3 < std::fabs(j1 - j2))
  {
    //G4cout << "j3 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m1) > j1)
  {
    //G4cout << "m1 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m2) > j2)
  {
    //G4cout << "m2 is out of bounds." << G4endl;
    return 0;
  }
  if(std::fabs(m3) > j3)
  {
    //G4cout << "m3 is out of bounds." << G4endl;
    return 0;
  }

  G4double out = (std::pow((-1),(j1-j2-m3)))/(std::pow((2*j3+1),(1.0/2.0)))*ClebschGordan(j1,m1,j2,m2,j3,-1*m3);
  return out;
}

// Rishita -- modified to include b parameters --> set to 1 for a gamma transition
G4double G4NuclearLevel::secondTransition(G4int k, G4double ji, G4double jf, G4int L1, G4int L2, G4double delta, G4double b1, G4double b2, G4double b3)
{
  G4double out;
  out = (1/(1+std::pow(delta,2)))*(b1*F(k,ji,L1,L1,jf)+2*delta*b2*F(k,ji,L1,L2,jf)+delta*delta*b3*F(k,ji,L2,L2,jf));
  return out;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
