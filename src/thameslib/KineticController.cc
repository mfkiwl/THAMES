/**
@file  KineticController.cc
@brief Method definitions for the KineticController class.

*/
#include "KineticController.h"

KineticController::KineticController() {
  temperature_ = 293.15;

  // default temperature (K)
  refT_ = 293.15;

  ///
  /// Clear out the vectors so they can be populated with values from the
  ///

  numPhases_ = 0;
  chemSys_ = NULL;
  lattice_ = NULL;
  phaseKineticModel_.clear();
  name_.clear();
  initScaledMass_.clear();
  scaledMass_.clear();
  specificSurfaceArea_.clear();
  refSpecificSurfaceArea_.clear();
  isKinetic_.clear();
  // waterId_ = 1;
  ICNum_ = 0;
  ICName_.clear();
  DCNum_ = 0;
  DCName_.clear();
  GEMPhaseNum_ = 0;

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion days or 27 million years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  verbose_ = warning_ = false;

  return;
}

KineticController::KineticController(ChemicalSystem *cs, Lattice *lattice,
                                     const string &fileName, const bool verbose,
                                     const bool warning)
    : chemSys_(cs), lattice_(lattice) {
  ///
  /// Clear out the vectors so they can be populated with values from the
  ///

  numPhases_ = 0;
  phaseKineticModel_.clear();
  name_.clear();
  isKinetic_.clear();

  // Set the verbose and warning flags

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
  cout << "KineticController::KineticController Constructor" << endl;
  cout.flush();
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  temperature_ = 293.15;
  refT_ = 293.15;

  ///
  /// Clear out the vectors so they can be populated with values from the
  /// JSON input file
  ///

  name_.clear();
  microPhaseId_.clear();

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion days or 27 million years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  ///
  /// Open the input JSON file for kinetic data and parse it
  ///

  string jsonext = ".json";
  size_t foundjson;
  foundjson = fileName.find(jsonext);
  try {
    if (foundjson != string::npos) {
      if (verbose_) {
        cout << "KineticModel data file is a JSON file" << endl;
      }
      parseDoc(fileName);
    } else {
      throw FileException("KineticModel", "KineticModel", fileName,
                          "NOT in JSON format");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  } catch (DataException dex) {
    dex.printException();
    exit(1);
  }

  int microPhaseId;

  if (verbose_) {
    cout << "KineticController::KineticController Finished reading "
            "chemistry.json "
         << endl;
    int size = microPhaseId_.size();
    for (int i = 0; i < size; ++i) {
      microPhaseId = microPhaseId_[i];
      if (isKinetic_[i]) {
        cout << "KineticController::KineticController kinetic phase "
             << microPhaseId << endl;
        cout << "KineticController::KineticController     name = "
             << chemSys_->getMicroPhaseName(microPhaseId) << endl;
      }
    }
    cout.flush();
  }

  // Assign the DC index for water

  // waterId_ = chemSys_->getDCId(WaterDCName);
  ICNum_ = chemSys_->getNumICs();
  DCNum_ = chemSys_->getNumDCs();
  ICName_ = chemSys_->getICName();
  DCName_ = chemSys_->getDCName();
  GEMPhaseNum_ = chemSys_->getNumGEMPhases();

  ICMoles_.resize(ICNum_, 0.0);
  ICMolesTot_.resize(ICNum_, 0.0);
  DCMoles_.resize(DCNum_, 0.0);
  DCMolesIni_.resize(DCNum_, 0.0);

  calcPhaseMasses();

  initScaledCementMass_ = chemSys_->getInitScaledCementMass();

  pKMsize_ = phaseKineticModel_.size();
  impurity_K2O_.resize(pKMsize_, 0);
  impurity_Na2O_.resize(pKMsize_, 0);
  impurity_Per_.resize(pKMsize_, 0);
  impurity_SO3_.resize(pKMsize_, 0);

  impurityDCID_.clear();
  impurityDCID_.push_back(chemSys_->getDCId("K2O"));
  impurityDCID_.push_back(chemSys_->getDCId("Na2O"));
  impurityDCID_.push_back(chemSys_->getDCId("Per")); // 170
  impurityDCID_.push_back(chemSys_->getDCId("SO3"));

  // initScaledMass_, scaledMass_ & scaledMassIni_ are
  // initialized in KineticController::parseMicroPhase :
  // initScaledMass_.push_back(0.0);
  // scaledMass_.push_back(0.0);
  // scaledMassIni_.push_back(0.0);

  return;
}

KineticController::~KineticController() {
  for (int midx = 0; midx < pKMsize_; ++midx) {
    delete phaseKineticModel_[midx];
  }
}

void KineticController::parseDoc(const string &docName) {
  int numEntry = -1; // Tracks number of solid phases
  int testgemid;

  ///
  /// The kineticData structure is used to temporarily hold parsed data
  /// for a given phase before the data are loaded permanently into class
  /// members.
  ///

  struct KineticData kineticData;

  /// Test for JSON existence

  ifstream f(docName.c_str());

  /// Parse the JSON file all at once
  json data = json::parse(f);
  f.close();

  try {

    /// Get an iterator to the root node of the JSON file
    /// @todo Add a better JSON validity check.

    json::iterator it = data.find("chemistry_data");
    json::iterator cdi = it.value().begin();

    // Test for non-emptiness
    if (cdi == it.value().end() || it == data.end()) {
      throw FileException("Controller", "parseDoc", docName, "Empty JSON file");
    }

    cdi = it.value().find("temperature");
    temperature_ = cdi.value();
    cdi = it.value().find("reftemperature");
    refT_ = cdi.value();
    cdi = it.value().find("phases");

    /// Each phase is a more complicated grouping of data that
    /// has a separate method for parsing.

    try {
      parseMicroPhases(cdi, numEntry, kineticData);
    } catch (DataException dex) {
      throw dex;
    }

    /// Push a copy of the isKinetic vector to the ChemicalSystem

    chemSys_->setIsKinetic(isKinetic_);

  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  /// All kinetic components have been parsed now.  Next, this block tries
  /// to handle pozzolanic effects (loi, SiO2 content, etc.) on any other
  /// kinetic phases

  setPozzEffectOnPK();

  return;
}

void KineticController::parseMicroPhases(const json::iterator cdi,
                                         int &numEntry,
                                         struct KineticData &kineticData) {
  int proposedgemphaseid, proposedDCid;
  int testgemid, testdcid;
  string testname;
  bool kineticfound = false;
  bool ispozz = false;
  bool isParrotKilloh = false;
  bool istherm = false;
  bool issol = false;

  json::iterator p = cdi.value().begin();

  for (int i = 0; i < cdi.value().size(); ++i) {
    initKineticData(kineticData);
    isKinetic_.push_back(false);
    p = cdi.value()[i].find("thamesname");
    testname = p.value();
    kineticData.name = testname;
    kineticData.microPhaseId = chemSys_->getMicroPhaseId(testname);
    kineticfound = ispozz = isParrotKilloh = false;

    p = cdi.value()[i].find("kinetic_data");
    if (p != cdi.value()[i].end()) {
      numEntry += 1;
      kineticfound = true;
      isKinetic_[isKinetic_.size() - 1] = true;
      kineticData.GEMPhaseId =
          chemSys_->getMicroPhaseToGEMPhase(kineticData.microPhaseId, 0);
      kineticData.DCId =
          chemSys_->getMicroPhaseDCMembers(kineticData.microPhaseId, 0);

      ///
      /// Kinetic data are grouped together,
      /// so there is a method written just for parsing that grouping
      ///

      try {
        parseKineticData(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    }
    if (kineticfound) {
      kineticData.scaledMass =
          chemSys_->getMicroPhaseMass(kineticData.microPhaseId);
      kineticData.temperature = temperature_;
      kineticData.reftemperature = refT_;
      makeModel(kineticData);
    }

    /// Some items should be added to vectors whether kinetically controlled or
    /// not

    name_.push_back(kineticData.name);
    microPhaseId_.push_back(kineticData.microPhaseId);
    initScaledMass_.push_back(0.0);
    scaledMass_.push_back(0.0);
    scaledMassIni_.push_back(0.0);
  }

  return;
}

void KineticController::parseKineticData(const json::iterator p,
                                         struct KineticData &kineticData) {
  bool typefound = false;

  try {
    json::iterator pp = p.value().find("type");
    kineticData.type = pp.value();
    if (kineticData.type == ParrotKillohType) {
      typefound = true;
      try {
        parseKineticDataForParrotKilloh(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    } else if (kineticData.type == StandardType) {
      typefound = true;
      try {
        parseKineticDataForStandard(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    } else if (kineticData.type == PozzolanicType) {
      typefound = true;
      try {
        parseKineticDataForPozzolanic(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    } else {
      throw HandleException("KineticController", "parseKineticData", "type",
                            "Model type not found");
    }

    if (!typefound) {
      throw HandleException("KineticController", "parseKineticData", "type",
                            "Model type not specified");
    }
  } catch (HandleException hex) {
    hex.printException();
  }

  return;
}

void KineticController::parseKineticDataForParrotKilloh(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    cout << "--->Parsing PK data for " << kineticData.name << endl;
    cout.flush();
  }

  json::iterator pp = p.value().find("specificSurfaceArea");
  kineticData.specificSurfaceArea = pp.value();

  // Reference specific surface area (m2/kg)
  pp = p.value().find("refSpecificSurfaceArea");
  kineticData.refSpecificSurfaceArea = pp.value();

  // Parrot-Killoh k1 parameter
  pp = p.value().find("k1");
  kineticData.k1 = pp.value();

  // Parrot-Killoh k2 parameter
  pp = p.value().find("k2");
  kineticData.k2 = pp.value();

  // Parrot-Killoh k3 parameter
  pp = p.value().find("k3");
  kineticData.k3 = pp.value();

  // Parrot-Killoh n1 parameter
  pp = p.value().find("n1");
  kineticData.n1 = pp.value();

  // Parrot-Killoh n3 parameter
  pp = p.value().find("n3");
  kineticData.n3 = pp.value();

  // Parrot-Killoh DOR_Hcoeff parameter
  pp = p.value().find("dorHcoeff");
  kineticData.dorHcoeff = pp.value();

  // Activation energy
  pp = p.value().find("activationEnergy");
  kineticData.activationEnergy = pp.value();

  return;
}

void KineticController::parseKineticDataForStandard(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    cout << "--->Parsing standard kinetic data for " << kineticData.name
         << endl;
    cout.flush();
  }

  // Specific surface area (m2/kg)
  json::iterator pp = p.value().find("specificSurfaceArea");
  if (pp != p.value().end()) {
    kineticData.specificSurfaceArea = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "specificSurfaceArea not found");
  }

  // Reference specific surface area (m2/kg)
  pp = p.value().find("refSpecificSurfaceArea");
  if (pp != p.value().end()) {
    kineticData.refSpecificSurfaceArea = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "refSpecificSurfaceArea not found");
  }

  // Dissolution rate constant (mol/m2/s)
  pp = p.value().find("dissolutionRateConst");
  if (pp != p.value().end()) {
    kineticData.dissolutionRateConst = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "dissolutionRateConst not found");
  }

  // Rate constant for early-age diffusion (mol/m2/s)
  pp = p.value().find("diffusionRateConstEarly");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstEarly = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "diffusionRateConstEarly not found");
  }

  // Dissolution rate constant (mol/m2/s)
  pp = p.value().find("diffusionRateConstLate");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstLate = pp.value();
  } else {
    kineticData.diffusionRateConstLate = kineticData.diffusionRateConstEarly;
    cout << "WARNING: For " << kineticData.name
         << " diffusionRateConstLate not found; setting it to "
            "diffusionRateConstEarly"
         << endl;
  }

  // Number of DC units produced in dissociation reaction
  pp = p.value().find("dissolvedUnits");
  if (pp != p.value().end()) {
    kineticData.dissolvedUnits = pp.value();
  } else {
    kineticData.dissolvedUnits = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dissolvedUnits not found; setting it to 1" << endl;
  }

  // Exponent on  the saturation index in the rate equation
  pp = p.value().find("siexp");
  if (pp != p.value().end()) {
    kineticData.siexp = pp.value();
  } else {
    kineticData.sio2 = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " sio2 not found; setting it to 1" << endl;
  }

  // Exponent on  the driving force term in the rate equation
  pp = p.value().find("dfexp");
  if (pp != p.value().end()) {
    kineticData.dfexp = pp.value();
  } else {
    kineticData.dfexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dfexp not found; setting it to 1" << endl;
  }

  // Loss on ignition of the material
  pp = p.value().find("loi");
  if (pp != p.value().end()) {
    kineticData.loi = pp.value();
  } else {
    kineticData.loi = 1.0e-6;
    cout << "WARNING: For " << kineticData.name
         << " loi not found; setting it to 1.0e-6" << endl;
  }

  // Activation energy for dissolution
  pp = p.value().find("activationEnergy");
  if (pp != p.value().end()) {
    kineticData.activationEnergy = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "activationEnergy not found");
  }

  return;
}

void KineticController::parseKineticDataForPozzolanic(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    cout << "--->Parsing pozzolanic data for " << kineticData.name << endl;
    cout.flush();
  }

  json::iterator pp = p.value().find("specificSurfaceArea");
  if (pp != p.value().end()) {
    kineticData.specificSurfaceArea = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "specificSurfaceArea not found");
  }

  pp = p.value().find("refSpecificSurfaceArea");
  if (pp != p.value().end()) {
    kineticData.refSpecificSurfaceArea = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "refSpecificSurfaceArea not found");
  }

  // Dissolution rate constant (mol/m2/s)
  pp = p.value().find("dissolutionRateConst");
  if (pp != p.value().end()) {
    kineticData.dissolutionRateConst = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "dissolutionRateConst not found");
  }

  // Early-age diffusion rate constant (mol/m2/s)
  pp = p.value().find("diffusionRateConstEarly");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstEarly = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "diffusionRateConstEarly not found");
  }

  // Later-age diffusion rate constant (mol/m2/s)
  pp = p.value().find("diffusionRateConstLate");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstLate = pp.value();
  } else {
    kineticData.diffusionRateConstLate = kineticData.diffusionRateConstEarly;
    cout << "WARNING: For " << kineticData.name
         << " diffusionRateConstLate not found; setting it to "
            "diffusionRateConstEarly"
         << endl;
  }

  // Number of DC units produced in dissociation reaction
  pp = p.value().find("dissolvedUnits");
  if (pp != p.value().end()) {
    kineticData.dissolvedUnits = pp.value();
  } else {
    kineticData.dissolvedUnits = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dissolvedUnits not found; setting it to 1" << endl;
  }

  // Exponent on the saturation index in the rate equation
  pp = p.value().find("siexp");
  if (pp != p.value().end()) {
    kineticData.siexp = pp.value();
  } else {
    kineticData.siexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " siexp not found; setting it to 1" << endl;
  }

  // Exponent on the driving force term in the rate equation
  pp = p.value().find("dfexp");
  if (pp != p.value().end()) {
    kineticData.dfexp = pp.value();
  } else {
    kineticData.dfexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dfexp not found; setting it to 1" << endl;
  }

  // Exponent on the degree of reaction term in the diffusion rate equation
  pp = p.value().find("dorexp");
  if (pp != p.value().end()) {
    kineticData.dorexp = pp.value();
  } else {
    kineticData.dorexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dorexp not found; setting it to 1" << endl;
  }

  // Exponent on the hydroxy ion activity in the rate equation
  pp = p.value().find("ohexp");
  if (pp != p.value().end()) {
    kineticData.ohexp = pp.value();
  } else {
    kineticData.ohexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " ohexp not found; setting it to 1" << endl;
  }

  // SiO2 mass fraction in the material
  pp = p.value().find("sio2");
  if (pp != p.value().end()) {
    kineticData.sio2 = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "sio2 not found");
  }

  // Al2O3 mass fraction in the material
  pp = p.value().find("al2o3");
  if (pp != p.value().end()) {
    kineticData.al2o3 = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "al2o3 not found");
  }

  // CaO mass fraction in the material
  pp = p.value().find("cao");
  if (pp != p.value().end()) {
    kineticData.cao = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "cao not found");
  }

  // Loss on ignition of the material
  pp = p.value().find("loi");
  if (pp != p.value().end()) {
    kineticData.loi = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "loi not found");
  }

  pp = p.value().find("activationEnergy");
  if (pp != p.value().end()) {
    kineticData.activationEnergy = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "activationEnergy not found");
  }

  return;
}

void KineticController::calcPhaseMasses(void) {
  int microPhaseId;
  double pscaledMass = 0.0;

  int size = microPhaseId_.size();

  for (int i = 0; i < size; i++) {
    microPhaseId = microPhaseId_[i];
    if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
      pscaledMass = chemSys_->getMicroPhaseMass(microPhaseId);
      scaledMass_[i] = pscaledMass;
      initScaledMass_[i] = pscaledMass;
      scaledMassIni_[i] = pscaledMass;

      // Setting the phase mass will also automatically calculate the phase
      // volume

      if (verbose_) {
        cout
            << "KineticController::getPhaseMasses reads solid micphase mass of "
            << chemSys_->getMicroPhaseName(microPhaseId) << " as "
            << initScaledMass_[i] << endl;
        cout.flush();
      }
    }
  }

  return;
}

double KineticController::getSolidMass(void) {
  int microPhaseId;
  double totmass = 0.0;
  int size = microPhaseId_.size();
  for (int i = 0; i < size; i++) {
    microPhaseId = microPhaseId_[i];
    if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
      totmass += chemSys_->getMicroPhaseMass(microPhaseId);
    }
  }

  return (totmass);
}

void KineticController::makeModel(struct KineticData &kineticData) {
  KineticModel *km = NULL;

  if (kineticData.type == ParrotKillohType) {
    // Read remaining Parrot and Killoh model parameters
    km = new ParrotKillohModel(chemSys_, lattice_, kineticData, verbose_,
                               warning_);
  } else if (kineticData.type == StandardType) {
    // Read remaining pozzolanic model parameters
    km = new StandardKineticModel(chemSys_, lattice_, kineticData, verbose_,
                                  warning_);
  } else if (kineticData.type == PozzolanicType) {
    // Read remaining pozzolanic model parameters
    km = new PozzolanicModel(chemSys_, lattice_, kineticData, verbose_,
                             warning_);
  }

  phaseKineticModel_.push_back(km);

  return;
}

void KineticController::setPozzEffectOnPK(void) {

  /// @todo This is the block where the influence of some components on the
  /// kinetic parameters of other components can be set.

  double refloi = 0.8;
  double loi = refloi;
  double maxloi = refloi;
  double fillareaeff = 1.0;
  double sio2val = 0.94;
  double refsio2val = 0.94;
  double betval = 29.0;
  double refbetval = 29.0;
  // double minpozzeffect = 1000.0;
  double minpozzeffect = 1.0;
  double pozzeffect = 1.0;

  int size = phaseKineticModel_.size();

  for (int midx = 0; midx < size; ++midx) {
    loi = phaseKineticModel_[midx]->getLossOnIgnition();
    if (loi > maxloi)
      maxloi = loi;
    if (phaseKineticModel_[midx]->getType() == PozzolanicType) {
      sio2val = phaseKineticModel_[midx]->getSio2();
      betval = phaseKineticModel_[midx]->getSpecificSurfaceArea();
      refbetval = phaseKineticModel_[midx]->getRefSpecificSurfaceArea();
      pozzeffect = pow((sio2val / refsio2val), 2.0) * (betval / refbetval);
      if (pozzeffect < minpozzeffect)
        minpozzeffect = pozzeffect;
      cout << endl
           << "KineticController::setPozzEffectOnPK for midx = " << midx
           << " (microPhaseId =  "
           << phaseKineticModel_[midx]->getMicroPhaseId()
           << ", microPhaseName = " << phaseKineticModel_[midx]->getName()
           << endl;

      cout << "  Ref LOI = " << refloi << endl;
      cout << "  LOI     = " << loi << endl;
      cout << "  Max LOI = " << maxloi << endl;
      cout << "  SiO2     = " << sio2val << endl;
      cout << "  Ref SiO2 = " << refsio2val << endl;
      cout << "  BET      = " << betval << endl;
      cout << "  Ref BET  = " << refbetval << endl;
      cout << "  Pozz Effect     = " << pozzeffect << endl;
      cout << "  Min Pozz Effect = " << minpozzeffect << endl;
      cout.flush();
    }
  }

  minpozzeffect *= (refloi / maxloi);

  /// The way this is set up, 0.0 <= refloi / maxloi <= 1.0
  for (int midx = 0; midx < size; ++midx) {
    if (phaseKineticModel_[midx]->getType() == ParrotKillohType) {
      phaseKineticModel_[midx]->setPfk(minpozzeffect);
    }
  }

  return;
}

void KineticController::calculateKineticStep(const double timestep, int cyc) {
  ///
  /// Initialize local variables
  ///
  ///

  int i;

  // double massDissolved = 0.0;
  cout << scientific << setprecision(15);
  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  // vector<int> impurityDCID;
  // impurityDCID.clear();
  // impurityDCID.push_back(chemSys_->getDCId("K2O"));
  // impurityDCID.push_back(chemSys_->getDCId("Na2O"));
  // impurityDCID.push_back(chemSys_->getDCId("Per")); // 170
  // impurityDCID.push_back(chemSys_->getDCId("SO3"));

  // cout << endl << "impurityDCID : " << endl;
  // for(i = 0; i < chemSys_->getNumMicroImpurities(); i++){
  //     cout << i << "\t" << impurityDCID[i] << endl; cout.flush();
  // }
  // cout << endl ;

  double totMassImpurity, massImpurity;

  int DCId;
  // int pKMsize = phaseKineticModel_.size();
  // static vector<double> scaledMassIni;
  double keepNumDCMoles;
  vector<int> phaseDissolvedId;
  phaseDissolvedId.resize(pKMsize_, 0);
  double numDCMolesDissolved, scaledMass, massDissolved;

  // static double hyd_time_ini = 0.0;
  double hyd_time = 0.0;

  for (i = 0; i < ICNum_; i++) {
    ICMoles_[i] = 0.0;
  }

  bool doTweak = (chemSys_->getTimesGEMFailed() > 0) ? true : false;

  if (doTweak) {
    hyd_time = hyd_time_ini_ + timestep;
    if (verbose_) {
      cout << endl
           << "  KineticController::calculateKineticStep - tweak before for "
              "cyc = "
           << cyc << "   hyd_time: " << hyd_time
           << "   hyd_time_ini_: " << hyd_time_ini_
           << "   timestep: " << timestep << endl;
    }
    for (int midx = 0; midx < pKMsize_; ++midx) {
      phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
      chemSys_->setMicroPhaseMass(phaseDissolvedId[midx], scaledMassIni_[midx]);
      if (verbose_) {
        cout << "    midx = " << midx
             << "     scaledMassIni[midx] = " << scaledMassIni_[midx] << endl;
      }
    }

    for (i = 0; i < DCNum_; i++) {
      DCMoles_[i] = DCMolesIni_[i];
    }

  } else {

    hyd_time = hyd_time_ini_ + timestep;
    cout << endl
         << "  KineticController::calculateKineticStep - initial for cyc = "
         << cyc << "   hyd_time: " << hyd_time
         << "   hyd_time_ini_: " << hyd_time_ini_ << "   timestep: " << timestep
         << endl;

    for (int midx = 0; midx < pKMsize_; ++midx) {
      phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
      scaledMassIni_[midx] =
          chemSys_->getMicroPhaseMass(phaseDissolvedId[midx]);
      if (verbose_) {
        cout << "    midx = " << midx
             << "     scaledMassIni[midx] = " << scaledMassIni_[midx] << endl;
      }
    }

    for (i = 0; i < DCNum_; i++) {
      DCMoles_[i] = chemSys_->getDCMoles(i);
      DCMolesIni_[i] = DCMoles_[i];
    }
  }

  if (verbose_) {
    cout << "  KineticController::calculateKineticStep Hydration Time = "
         << hyd_time << endl;
    cout.flush();
  }

  try {
    // cout << "  KineticController::calculateKineticStep     hyd_time = "
    //      << hyd_time << "\tcyc = " << cyc << endl;

    if (hyd_time < leachTime_ && hyd_time < sulfateAttackTime_) {

      // if (!doTweak) {
      //  @todo BULLARD PLACEHOLDER
      //  Still need to implement constant gas phase composition
      //  Will involve equilibrating gas with aqueous solution
      //
      //  First step each iteration is to equilibrate gas phase
      //  with the electrolyte, while forbidding anything new
      //  from precipitating.

      /// This is a big kluge for internal relative humidity
      /// @note Using new gel and interhydrate pore size distribution model
      ///       which is currently contained in the Lattice object.
      ///
      /// Surface tension of water is gamma = 0.072 J/m2
      /// Molar volume of water is Vm = 1.8e-5 m3/mole
      /// The Kelvin equation is
      ///    p/p0 = exp (-4 gamma Vm / d R T) = exp (-6.23527e-7 / (d T))
      ///
      ///    where d is the pore diameter in meters and T is absolute
      ///    temperature

      /// Assume a zero contact angle for now.
      /// @todo revisit the contact angle issue

      /// Loop over all kinetic models

      //*******
      double totalDOR = 0;

      if (initScaledCementMass_ > 0)
        totalDOR = (initScaledCementMass_ - chemSys_->getScaledCementMass()) /
                   initScaledCementMass_;

      double dcmoles;
      for (int midx = 0; midx < pKMsize_; ++midx) {

        DCId = phaseKineticModel_[midx]->getDCId();
        phaseKineticModel_[midx]->calculateKineticStep(
            timestep, scaledMass, massDissolved, cyc, totalDOR);

        // chemSys_->setMicroPhaseMass(phaseDissolvedId[midx], scaledMass);
        // chemSys_->setMicroPhaseMassDissolved(phaseDissolvedId[midx],
        //                                      massDissolved);
        chemSys_->updateMicroPhaseMasses(phaseDissolvedId[midx], scaledMass, 0);

        if (verbose_) {
          cout << "New scaled mass = "
               << chemSys_->getMicroPhaseMass(phaseDissolvedId[midx])
               << " and new volume = "
               << chemSys_->getMicroPhaseVolume(phaseDissolvedId[midx]) << endl;
          cout.flush();
        }

        totMassImpurity = 0;
        keepNumDCMoles = 0;
        numDCMolesDissolved = 0;

        massImpurity = massDissolved * chemSys_->getK2o(phaseDissolvedId[midx]);
        totMassImpurity += massImpurity;
        dcmoles = massImpurity / chemSys_->getDCMolarMass("K2O");
        DCMoles_[impurityDCID_[0]] += dcmoles;
        impurity_K2O_[midx] = dcmoles;
        // DCMoles_[impurityDCID[0]] +=
        //     massImpurity / chemSys_->getDCMolarMass("K2O");

        massImpurity =
            massDissolved * chemSys_->getNa2o(phaseDissolvedId[midx]);
        totMassImpurity += massImpurity;
        dcmoles = massImpurity / chemSys_->getDCMolarMass("Na2O");
        DCMoles_[impurityDCID_[1]] += dcmoles;
        impurity_Na2O_[midx] = dcmoles;
        // DCMoles_[impurityDCID[1]] +=
        //     massImpurity / chemSys_->getDCMolarMass("Na2O");

        massImpurity = massDissolved * chemSys_->getMgo(phaseDissolvedId[midx]);
        totMassImpurity += massImpurity;
        dcmoles = massImpurity / chemSys_->getDCMolarMass("Per");
        DCMoles_[impurityDCID_[2]] += dcmoles;
        impurity_Per_[midx] = dcmoles;
        // DCMoles_[impurityDCID[2]] +=
        //     massImpurity / chemSys_->getDCMolarMass("Per"); // MgO

        massImpurity = massDissolved * chemSys_->getSo3(phaseDissolvedId[midx]);
        totMassImpurity += massImpurity;
        dcmoles = massImpurity / chemSys_->getDCMolarMass("SO3");
        DCMoles_[impurityDCID_[3]] += dcmoles;
        impurity_SO3_[midx] = dcmoles;
        // DCMoles_[impurityDCID[3]] +=
        //     massImpurity / chemSys_->getDCMolarMass("SO3");

        numDCMolesDissolved =
            (massDissolved - totMassImpurity) / chemSys_->getDCMolarMass(DCId);
        keepNumDCMoles = DCMoles_[DCId] - numDCMolesDissolved;
        chemSys_->setDCLowerLimit(DCId, keepNumDCMoles);
        if (verbose_) {
          cout << "    calculateKineticStep - "
                  "midx/DCId/DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
               << midx << " / " << DCId << " / " << DCMoles_[DCId] << " / "
               << numDCMolesDissolved << " / " << keepNumDCMoles << endl;
          cout << "    calculateKineticStep - scaledMass/massDissolved/"
                  "totMassImpurity/massDissolved - totMassImpurity : "
               << scaledMass << " / " << massDissolved << " / "
               << totMassImpurity << " / " << massDissolved - totMassImpurity
               << endl;
        }
      }

      if (verbose_ && doTweak) {
        cout << endl
             << "  KineticController::calculateKineticStep "
                "- tweak after for cyc = "
             << cyc << endl;
      }

    } // End of normal hydration block
  } // End of try block

  catch (EOBException eex) {
    eex.printException();
    exit(1);
  } catch (DataException dex) {
    dex.printException();
    exit(1);
  } catch (FloatException fex) {
    fex.printException();
    exit(1);
  } catch (out_of_range &oor) {
    EOBException ex("KineticController", "calculateKineticStep", oor.what(), 0,
                    0);
    ex.printException();
    exit(1);
  }

  for (i = 0; i < DCNum_; i++) {
    // cout << " " << i << "\t" << DCName_[i] << ": " << DCMoles_[i] << " mol"
    // << endl;
    chemSys_->setDCMoles(i, DCMoles_[i]);
    // cout << "          " << DCName_[i] << ": " << chemSys_->getDCMoles(i) <<
    // " mol" << endl;
  }

  if (verbose_) {
    cout << "  KineticController::calculateKineticStep end - cyc = " << cyc
         << endl;
    cout.flush();
  }

  return;
}

double KineticController::updateKineticStep(int cyc, int pId,
                                            double scaledMass) {
  int i;
  string modelName;
  double totMassImpurity, massImpurity;
  int DCId;
  double keepNumDCMoles;
  int phaseDissolvedId;
  double numDCMolesDissolved, massDissolved;

  for (i = 0; i < ICNum_; i++) {
    ICMoles_[i] = 0.0;
  }

  int midx;
  for (midx = 0; midx < pKMsize_; ++midx) {
    phaseDissolvedId = phaseKineticModel_[midx]->getMicroPhaseId();
    if (pId == phaseDissolvedId) {
      DCId = phaseKineticModel_[midx]->getDCId();
      break;
    }
  }
  chemSys_->setMicroPhaseMass(phaseDissolvedId, scaledMassIni_[midx]);
  modelName =
      phaseKineticModel_[midx]
          ->getModelName(); // updateKineticStep(scaledMass , massDissolved);
  cout << endl
       << "  KineticController::updateKineticStep - for cyc = " << cyc
       << " & phaseId = " << pId << endl;
  cout << "    midx = " << midx << "   modelName : " << modelName
       << "   scaledMassIni[midx] = " << scaledMassIni_[midx]
       << "   scaledMass = " << scaledMass << endl;

  // int numMicroPhases = isKinetic_.size();
  // int microPhaseMemberSize;
  // vector<int> DCIdMicroPhaseMembers;
  // for (i = 1; i < numMicroPhases; i++) {
  //   if (!isKinetic_[i]) {
  //     DCIdMicroPhaseMembers = chemSys_->getMicroPhaseDCMembers(i);
  //     microPhaseMemberSize = DCIdMicroPhaseMembers.size();
  //     for (int j = 0; j < microPhaseMemberSize; j++) {
  //        DCMoles_[DCIdMicroPhaseMembers[j]] =
  //        DCMolesIni_[DCIdMicroPhaseMembers[j]];
  //     }
  //   }
  // }

  for (i = 0; i < DCNum_; i++) {
    if (chemSys_->getIsDCKinetic(i)) {
      if (i == DCId)
        DCMoles_[DCId] = DCMolesIni_[DCId];
    } else {
      if ((i != impurityDCID_[0]) && (i != impurityDCID_[1]) &&
          (i != impurityDCID_[2]) && (i != impurityDCID_[3]))
        DCMoles_[i] = DCMolesIni_[i];
    }
  }
  // DCMoles_[DCId] = DCMolesIni_[DCId];
  DCMoles_[impurityDCID_[0]] -= impurity_K2O_[midx];
  DCMoles_[impurityDCID_[1]] -= impurity_Na2O_[midx];
  DCMoles_[impurityDCID_[2]] -= impurity_Per_[midx];
  DCMoles_[impurityDCID_[3]] -= impurity_SO3_[midx];

  /// for this kinetic model

  massDissolved = scaledMassIni_[midx] - scaledMass;

  // chemSys_->setMicroPhaseMass(phaseDissolvedId, scaledMass);
  // chemSys_->setMicroPhaseMassDissolved(phaseDissolvedId, massDissolved);
  chemSys_->updateMicroPhaseMasses(phaseDissolvedId, scaledMass, 1);

  totMassImpurity = 0;
  keepNumDCMoles = 0;
  numDCMolesDissolved = 0;

  double dcmoles;
  massImpurity = massDissolved * chemSys_->getK2o(phaseDissolvedId);
  totMassImpurity += massImpurity;
  dcmoles = massImpurity / chemSys_->getDCMolarMass("K2O");
  DCMoles_[impurityDCID_[0]] += dcmoles;
  impurity_K2O_[midx] = dcmoles;
  // DCMoles_[impurityDCID[0]] +=
  //     massImpurity / chemSys_->getDCMolarMass("K2O");

  massImpurity = massDissolved * chemSys_->getNa2o(phaseDissolvedId);
  totMassImpurity += massImpurity;
  dcmoles = massImpurity / chemSys_->getDCMolarMass("Na2O");
  DCMoles_[impurityDCID_[1]] += dcmoles;
  impurity_Na2O_[midx] = dcmoles;
  // DCMoles_[impurityDCID[1]] +=
  //     massImpurity / chemSys_->getDCMolarMass("Na2O");

  massImpurity = massDissolved * chemSys_->getMgo(phaseDissolvedId);
  totMassImpurity += massImpurity;
  dcmoles = massImpurity / chemSys_->getDCMolarMass("Per");
  DCMoles_[impurityDCID_[2]] += dcmoles;
  impurity_Per_[midx] = dcmoles;
  // DCMoles_[impurityDCID[2]] +=
  //     massImpurity / chemSys_->getDCMolarMass("Per"); // MgO

  massImpurity = massDissolved * chemSys_->getSo3(phaseDissolvedId);
  totMassImpurity += massImpurity;
  dcmoles = massImpurity / chemSys_->getDCMolarMass("SO3");
  DCMoles_[impurityDCID_[3]] += dcmoles;
  impurity_SO3_[midx] = dcmoles;
  // DCMoles_[impurityDCID[3]] +=
  //     massImpurity / chemSys_->getDCMolarMass("SO3");

  numDCMolesDissolved =
      (massDissolved - totMassImpurity) / chemSys_->getDCMolarMass(DCId);
  keepNumDCMoles = DCMoles_[DCId] - numDCMolesDissolved;
  chemSys_->setDCLowerLimit(DCId, keepNumDCMoles);
  cout << "    updateKineticStep "
          "midx/DCId/DCMoles_/numDCMolesDissolved/keepNumDCMoles[midx]/"
          "massDissolved/"
          "totMassImpurity/massDissolved - totMassImpurity : "
       << midx << " / " << DCId << " / " << DCMoles_[DCId] << " / "
       << numDCMolesDissolved << " / " << keepNumDCMoles << " / "
       << massDissolved << " / " << totMassImpurity << " / "
       << massDissolved - totMassImpurity << endl;

  for (i = 0; i < DCNum_; i++) {
    // cout << " " << i << "\t" << DCName_[i] << ": " << DCMoles_[i] << " mol"
    // << endl;
    chemSys_->setDCMoles(i, DCMoles_[i]);
    // cout << "          " << DCName_[i] << ": " << chemSys_->getDCMoles(i) <<
    // " mol" << endl;
  }

  cout << "  KineticController::updateKineticStep end - cyc = " << cyc << endl;
  cout.flush();

  return massDissolved;
}
