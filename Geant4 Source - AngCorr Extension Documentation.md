# Geant4 Gamma-Electron Angular Correlations Extension
### Documentation
----------
##### Modified Modules: December 2017
##### Rishita Gudapati  

* _"Evan's version" refers to the Geant4 Gamma-Gamma Angular Correlations extension written by Evan Rand.
      It is available on the GRIFFINCollaboration GitHub page under `Geant4GammaGammaAngularCorrelations10.01.p01`; this is         the version I started working from._
* _"original source files" refers to the Geant4 version 10.01.p03 source files._

### G4DiscreteGammaTransition
----
_This is where the program decides whether to emit a $$\gamma$$ or an $$e^-$$ based on the internal conversion coefficients in the Photon Evaporation file. The actual generation of the particle is done in `G4VGammaDeexcitation`._
 
 _This is also where `GetThetaFromWTheta` is located (called in `G4VGammaDeexcitation`). This method returns a $$\theta$$ value consistent with a specific correlation (passed in as casenumber). The correlation is not applied here._
 
 ##### Header File
 Additions to Evan's version:
 * modify `GetThetaFromWTheta` to accept casenumber argument
 * declarations for vectors mentioned below
 
 ##### Source File
 Additions to Evan's version:
 * each case (gg, ge..) has its own set of coefficients for which the vectors are populated at each level
 * the above is also true for the _maxWTheta_ that belongs to each case
 * remove the subtraction of bond energy for electrons: this has been moved to `G4VGammaDeexcitation`
		to be consistent with the newer _Geant4_ versions
 * normalization is now performed separately for each case
   * I'm still not entirely sure how this normalization procedure works, and it is the likely culprit behind our
	   inability to simulate the decay of _Tl198[543.6]_ -- see the emails from Evan, which Mike Bowry has copies of
 * the normalization procedure produces a theta value which is then returned to `G4VGammaDeexcitation` and applied
	  to the polarization vector
 


### G4NuclearLevel
----
_This is where the values of the angular coefficients are calculated._
_This is also where `GenerateWThetaParameters` is located._

##### Header File
Additions to Evan's version: 
* declarations for vectors mentioned below
* declare methods for filling _b_, $$\alpha$$ parameters for a specific level
	  used in `G4NuclearLevelManager` when reading and storing values in multipole file

##### Source File
Additions to Evan's version:
* create and initialize vectors for _b~1~_ through _b~6~_ (_b~1~_ to _b~3~_ belong to _a~2~_ coefficient, _b~4~_ to _b~6~_ belong to _a~4~_)
* create and initialize vectors for _a~2~_ to _a~10~_ for *each* separate transition type (gg, ge, eg, ee)
* create and initialize vectors for _maxWTheta_ for each separate transition type
* calculate mixing ratios for $$e^-$$ transitions
* coefficients: modify first and second transition equations to include _b_, $$\alpha$$ parameters
  * _b_ parameters are set to 1 for a $$\gamma$$ transition
* modified `GenerateWThetaParameters` to accept _b_, $$\alpha$$ values given in multipole file 
  * this function is called in `G4NuclearLevelManager` where the multipole file is read in
* calculate and output angular coefficients for all transition types
* suppress output for levels with _a~2~($$\gamma\gamma$$)_ $$=$$ _a~4~($$\gamma\gamma$$)_ $$= 0$$
  * we were using the entire Photon Evaporation file. In theory, this shouldn't be
	   necessary if we only include the levels that we're interested in
* calculate _maxWTheta_ for all four cases -- used for normalization purposes in `G4DiscreteGammaTransition`

### G4NuclearLevelManager
----
_This is where we read in the multipole file and save its data to the appropriate energy levels_
_This is also where `GenerateWThetaParameters` is called_
##### Header File
Additions to Evan's version:
* declarations for _b_ and $$\alpha$$ vectors 
  * each vector consists of one element for each energy level in the Photon Evaporation file

##### Source File
Additions to Evan's version:
* read in _b_ parameters and $$\alpha$$ values from multipole file
* fill _b_ parameters and $$\alpha$$ values into the right level 
  * it finds the right level by looping through all levels in the Photon Evaporation file and comparing the level energy in the multipole file to the current level energy
  * in checking that the level is compared correctly, I ran into some floating point precision errors, so I
	     added the _firstcheck_ and _secondcheck_ variables to counteract them. The precision can probably increased
	     further in these checks.
* retrieve the _b_ and $$\alpha$$ values for each level in the Photon Evaporation file (two levels at a time -- the 	  current one and the one immediately above it in energy) and pass these values to `G4NuclearLevel::GenerateWThetaParameters`

### G4NuclearLevelStore
----
_This is just a storage module for variables used across the code_

##### Header File
Deletions from Evan's version:
* removed getter and setter for _FirstGammaDecay_

Additions to Evan's version:
* getter and setter for _G4bool FirstProduct_
  * every nuclear deexcitation results in a product chain inside the `DoDecay` method of `G4VGammaDeexcitation`. 
  _FirstProduct_ simply stores whether or not the product list is empty. The first product of a deexcitation is isotropically emitted from the nucleus. 
* getter and setter for _G4bool ElecFirst_
  * true if the particle generated first in the cascade was an $$e^-$$
  * its value is only relevant for the second particle in the cascade
* getter and setter for _G4int CaseNumber_
  * check the current particle type and the previous particle type to determine case.
				(this is done in `G4VGammaDeexcitation`) 
  * $$>= 0$$ : only for second and subsequent particles
    * case 0: $$\gamma$$ -- $$\gamma$$
    * case 1: $$\gamma$$ -- $$e^-$$ (previous particle was a gamma, current particle is an electron)
    * case 2: $$e^-$$-- $$\gamma$$
    * case 3: $$e^-$$ -- $$e^-$$
  * $$-2$$ : for first particle (emitted isotropically)
  * $$-3$$ : for _NULL_ transitions
  * $$-1$$ : for transitions which are not $$\gamma$$ or $$e^-$$ AND are not the first product 
    * emitted isotropically 

##### Source File
Unchanged from Evan's version.

### G4VGammaDeexcitation
----
_This is where new particles in the cascade are generated_
_This is also where `GetThetaFromWTheta` is called and the correlation is applied to the newly generated particle_

##### Header File
Unchanged from Evan's version

##### Source File
Additions to Evan's version:
* `DoDecay` method: 
  * placed a check here to see if the product list is empty or not. This tells us whether  or not the particle being generated is the first in the cascade (stored as _G4bool FirstProduct_ in `G4NuclearLevelStore`)
* Case Selection
  * see _G4int CaseNumber_ above
  * determine which case the currently generated particle belongs to, and pass this number to 
	  `G4DiscreteGammaTransition::GetThetaFromWTheta`
  * first particle in cascade is emitted isotropically
  * all non-gamma and non-electron particles are emitted isotropically
* subtract bond energy of electrons here (INSTEAD of in `G4DiscreteGammaTransition`)
* takes the $$\theta$$ produced from `GetThetaFromWTheta` and applies it to a polarization vector to simulate the correlation

### G4RadioactiveDecay
----
_This is where radioactive decay modes are accessed and the correct decay path is chosen for the source particle_

The entire module is unchanged from the original _Geant4_ source files. Some deletions where made from Evan's version to get it back to this original state (removed _FirstGammaDecay_ checks). This is only included in the package in case the user has Evan's extension already installed.

### G4RadioactiveDecayMessenger
----
_Not entirely sure what this module does._

The entire module is unchanged from Evan's version. Only included in this package so that the user does not have to download both the Gamma-Gamma *and* the Gamma-Electron extensions.