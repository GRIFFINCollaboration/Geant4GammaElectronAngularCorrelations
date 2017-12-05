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
// $Id: G4NuclearLevel.hh 86986 2014-11-21 13:00:05Z gcosmo $
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
//      Creation date: 25 October 1998
//
//      Modifications:
//	  06 Oct 2010, M. Kelsey (kelsey@slac.stanford.edu)
//		Add friendship for G4NuclearLevelManager; define private
//		constructors without vectors.
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//
//        28 October 2010, V.Ivanchenko moved copy constructor to source, cleanup
//
//
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVEL_HH
#define G4NUCLEARLEVEL_HH 1

#include "globals.hh"
#include <vector>

class G4NuclearLevel
{

public:
  G4NuclearLevel(G4double energy, G4double halfLife,
         G4double angularMomentum, const std::vector<G4double>& eGamma,
         const std::vector<G4double>& wGamma, const std::vector<G4double>& polarities,
         const std::vector<G4double>& kCC, const std::vector<G4double>& l1CC,
         const std::vector<G4double>& l2CC, const std::vector<G4double>& l3CC,
         const std::vector<G4double>& m1CC, const std::vector<G4double>& m2CC,
         const std::vector<G4double>& m3CC, const std::vector<G4double>& m4CC,
         const std::vector<G4double>& m5CC, const std::vector<G4double>& nPlusCC,
         const std::vector<G4double>& totalCC);

  ~G4NuclearLevel();

  const std::vector<G4double>& GammaEnergies() const;

  const std::vector<G4double>& GammaWeights() const;

  const std::vector<G4double>& GammaProbabilities() const;

  const std::vector<G4double>& GammaCumulativeProbabilities() const;

  const std::vector<G4double>& GammaPolarities() const;

  const std::vector<G4double>& KConvertionProbabilities() const;

  const std::vector<G4double>& L1ConvertionProbabilities() const;

  const std::vector<G4double>& L2ConvertionProbabilities() const;

  const std::vector<G4double>& L3ConvertionProbabilities() const;

  const std::vector<G4double>& M1ConvertionProbabilities() const;

  const std::vector<G4double>& M2ConvertionProbabilities() const;

  const std::vector<G4double>& M3ConvertionProbabilities() const;

  const std::vector<G4double>& M4ConvertionProbabilities() const;

  const std::vector<G4double>& M5ConvertionProbabilities() const;

  const std::vector<G4double>& NPlusConvertionProbabilities() const;

  const std::vector<G4double>& TotalConvertionProbabilities() const;

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  const std::vector<G4double>& L1() const;

  const std::vector<G4double>& L2() const;

  const std::vector<G4double>& MixingRatio() const;
  
  const std::vector< std::vector< G4double > >& HigherLevelEnergy() const;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Rishita Gudapati- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  const std::vector< std::vector< G4double > >& A2_GG() const;

  const std::vector< std::vector< G4double > >& A4_GG() const;

  const std::vector< std::vector< G4double > >& A6_GG() const;

  const std::vector< std::vector< G4double > >& A8_GG() const;

  const std::vector< std::vector< G4double > >& A10_GG() const;

  const std::vector< std::vector< G4double > >& MaxWTheta_GG() const;

  const std::vector< std::vector< G4double > >& A2_GE() const;

  const std::vector< std::vector< G4double > >& A4_GE() const;

  const std::vector< std::vector< G4double > >& A6_GE() const;

  const std::vector< std::vector< G4double > >& A8_GE() const;

  const std::vector< std::vector< G4double > >& A10_GE() const;

  const std::vector< std::vector< G4double > >& MaxWTheta_GE() const;

  const std::vector< std::vector< G4double > >& A2_EG() const;

  const std::vector< std::vector< G4double > >& A4_EG() const;

  const std::vector< std::vector< G4double > >& A6_EG() const;

  const std::vector< std::vector< G4double > >& A8_EG() const;

  const std::vector< std::vector< G4double > >& A10_EG() const;

  const std::vector< std::vector< G4double > >& MaxWTheta_EG() const;

  const std::vector< std::vector< G4double > >& A2_EE() const;

  const std::vector< std::vector< G4double > >& A4_EE() const;

  const std::vector< std::vector< G4double > >& A6_EE() const;

  const std::vector< std::vector< G4double > >& A8_EE() const;

  const std::vector< std::vector< G4double > >& A10_EE() const;

  const std::vector< std::vector< G4double > >& MaxWTheta_EE() const;

  const std::vector< G4double >& B1() const;

  const std::vector< G4double >& B2() const;

  const std::vector< G4double >& B3() const;

  const std::vector< G4double >& B4() const;

  const std::vector< G4double >& B5() const;

  const std::vector< G4double >& B6() const;

  const std::vector< G4double >& AlphL() const;

  const std::vector< G4double >& AlphLp() const;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  G4double Energy() const;

  G4double AngularMomentum() const;

  G4double HalfLife() const;

  G4int NumberOfGammas() const;

  void PrintAll() const;

  void PrintLevels() const;

  G4bool operator==(const G4NuclearLevel &right) const;
  G4bool operator!=(const G4NuclearLevel &right) const;
  G4bool operator<(const G4NuclearLevel &right) const;
  G4NuclearLevel& operator=(const G4NuclearLevel &right);

  G4NuclearLevel(const G4NuclearLevel &right);

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  inline void FillL1(G4int i, G4double value) { _L1[i] = value; };
  inline void FillL2(G4int i, G4double value) { _L2[i] = value; };
  inline void FillMixingRatio(G4int i, G4double value) { _mixingRatio[i] = value; };
  
  void FillHigherLevelEnergy(G4int gamma_i, G4int hlevel_i, G4double hlevel_energy);
  
  // Rishita -- modified to include b parameters
  void GenerateWThetaParameters(G4int gamma_i, G4int hlevel_i, G4double hlevel_energy, G4double lowerGammaEnergy, G4double higherGammaEnergy, G4double ji, G4double jo, G4double jf, G4int L1, G4int L1p, G4int L2, G4int L2p, G4double delta1, G4double delta2, G4bool boolGoodLevelToOutputToScreen, G4double bf1, G4double bf2, G4double bf3, G4double bf4, G4double bf5, G4double bf6, G4double alphL1, G4double alphL1p, G4double bs1, G4double bs2, G4double bs3, G4double bs4, G4double bs5, G4double bs6, G4double alphL2, G4double alphL2p);

  G4double ClebschGordan(G4double j1, G4double m1, G4double j2, G4double m2, G4double j, G4double m);
  G4double Factorial(G4double value);
  G4double LegendreP(G4int n, G4double x);
  G4double F(G4int k, G4double jf, G4int L1, G4int L2, G4double ji);

  // Rishita -- modified to include b parameters
  G4double firstTransition (G4int k, G4double ji, G4double jf, G4int L1, G4int L2, G4double delta, G4double b1, G4double b2, G4double b3);

  G4double RacahW(G4double a, G4double b, G4double c, G4double d, G4double e, G4double f);
  G4double Wigner6j(G4double J1, G4double J2, G4double J3, G4double J4, G4double J5, G4double J6);
  G4double Wigner3j(G4double j1, G4double j2, G4double j3, G4double m1, G4double m2, G4double m3);

  // Rishita -- modified to include b parameters
  G4double secondTransition(G4int k, G4double ji, G4double jf, G4int L1, G4int L2, G4double delta, G4double b1, G4double b2, G4double b3);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  // Rishita Gudapati- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  inline void FillB1(G4int i, G4double value) {_b1[i] = value; };
  inline void FillB2(G4int i, G4double value) {_b2[i] = value; };
  inline void FillB3(G4int i, G4double value) {_b3[i] = value; };
  inline void FillB4(G4int i, G4double value) {_b4[i] = value; };
  inline void FillB5(G4int i, G4double value) {_b5[i] = value; };
  inline void FillB6(G4int i, G4double value) {_b6[i] = value; };
  inline void FillAlphL(G4int i, G4double value) {_alphL[i] = value; };
  inline void FillAlphLp(G4int i, G4double value) {_alphLp[i] = value; };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

private:

  friend class G4NuclearLevelManager;

  G4NuclearLevel();

  G4NuclearLevel(G4double energy, G4double halfLife,
         G4double angularMomentum);

  void Finalize();

  void MakeProbabilities();
  void MakeCumProb();
  void MakeMultipoleData(); // Evan Rand

  G4int Increment(G4int aF);

  std::vector<G4double> _energies;
  std::vector<G4double> _weights;
  std::vector<G4double> _prob;
  std::vector<G4double> _cumProb;
  std::vector<G4double> _polarities;
  std::vector<G4double> _kCC;
  std::vector<G4double> _l1CC;
  std::vector<G4double> _l2CC;
  std::vector<G4double> _l3CC;
  std::vector<G4double> _m1CC;
  std::vector<G4double> _m2CC;
  std::vector<G4double> _m3CC;
  std::vector<G4double> _m4CC;
  std::vector<G4double> _m5CC;
  std::vector<G4double> _nPlusCC;
  std::vector<G4double> _totalCC;

  G4double _energy;
  G4double _halfLife;
  G4double _angularMomentum;
  G4int _nGammas;

  // Evan Rand - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  std::vector<G4double> _L1;
  std::vector<G4double> _L2;
  std::vector<G4double> _mixingRatio;
  std::vector< std::vector< G4double > > _higherLevelEnergy;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  // Rishita Gudapati- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  std::vector<G4double> _b1;
  std::vector<G4double> _b2;
  std::vector<G4double> _b3;
  std::vector<G4double> _b4;
  std::vector<G4double> _b5;
  std::vector<G4double> _b6;
  std::vector<G4double> _alphL1;
  std::vector<G4double> _alphL2;

  std::vector< std::vector< G4double > > _a2_gg;
  std::vector< std::vector< G4double > > _a4_gg;
  std::vector< std::vector< G4double > > _a6_gg;
  std::vector< std::vector< G4double > > _a8_gg;
  std::vector< std::vector< G4double > > _a10_gg;
  std::vector< std::vector< G4double > > _maxWTheta_gg;

  std::vector< std::vector< G4double > > _a2_ge;
  std::vector< std::vector< G4double > > _a4_ge;
  std::vector< std::vector< G4double > > _a6_ge;
  std::vector< std::vector< G4double > > _a8_ge;
  std::vector< std::vector< G4double > > _a10_ge;
  std::vector< std::vector< G4double > > _maxWTheta_ge;

  std::vector< std::vector< G4double > > _a2_eg;
  std::vector< std::vector< G4double > > _a4_eg;
  std::vector< std::vector< G4double > > _a6_eg;
  std::vector< std::vector< G4double > > _a8_eg;
  std::vector< std::vector< G4double > > _a10_eg;
  std::vector< std::vector< G4double > > _maxWTheta_eg;

  std::vector< std::vector< G4double > > _a2_ee;
  std::vector< std::vector< G4double > > _a4_ee;
  std::vector< std::vector< G4double > > _a6_ee;
  std::vector< std::vector< G4double > > _a8_ee;
  std::vector< std::vector< G4double > > _a10_ee;
  std::vector< std::vector< G4double > > _maxWTheta_ee;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

};

#endif
