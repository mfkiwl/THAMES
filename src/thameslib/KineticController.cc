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
  waterId_ = 1;
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
  /// XML input file
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
  /// Open the input XML file for kinetic data and parse it
  ///

  string xmlext = ".xml";
  size_t foundxml;
  foundxml = fileName.find(xmlext);
  try {
    if (foundxml != string::npos) {
      if (verbose_) {
        cout << "KineticModel data file is an XML file" << endl;
      }
      parseDoc(fileName);
    } else {
      throw FileException("KineticModel", "KineticModel", fileName,
                          "NOT in XML format");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  int microPhaseId;

  if (verbose_) {
    cout << "KineticController::KineticController Finished reading "
            "chemistry.xml "
         << endl;
    for (int i = 0; i < microPhaseId_.size(); ++i) {
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

  waterId_ = chemSys_->getDCId(WaterDCName);
  ICNum_ = chemSys_->getNumICs();
  DCNum_ = chemSys_->getNumDCs();
  ICName_ = chemSys_->getICName();
  DCName_ = chemSys_->getDCName();
  GEMPhaseNum_ = chemSys_->getNumGEMPhases();

  ICMoles_.resize(ICNum_,0.0);
  ICMolesTot_.resize(ICNum_,0.0);
  DCMoles_.resize(DCNum_,0.0);

  calcPhaseMasses();

  initScaledCementMass_ = chemSys_->getInitScaledCementMass();

  return;
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

  ///
  /// This method uses the libxml library, so it needs to be added and linked
  /// at compile time.
  ///

  xmlDocPtr doc;
  xmlChar *key;
  xmlNodePtr cur;

  cout.flush();
  doc = xmlParseFile(docName.c_str());

  ///
  /// Check if the xml file is valid and parse it if so.
  /// @note This block requires schema file to be local

  try {
    string rxcsd = "chemistry.xsd";
    if (!is_xml_valid(doc, rxcsd.c_str())) {
      throw FileException("KineticModel", "KineticModel", docName,
                          "xml NOT VALID");
    }

    if (doc == NULL) {
      throw FileException("KineticModel", "KineticModel", docName,
                          "xml NOT parsed successfully");
    }

    cur = xmlDocGetRootElement(doc);

    if (cur == NULL) {
      xmlFreeDoc(doc);
      throw FileException("KineticModel", "KineticModel", docName,
                          "xml document is empty");
    }

    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"temperature"))) {
        key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        string st((char *)key);
        from_string(temperature_, st);
        xmlFree(key);
      } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"reftemperature"))) {
        key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        string st((char *)key);
        from_string(refT_, st);
        xmlFree(key);
      } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"phase"))) {

        /// Each phase is a more complicated grouping of data that
        /// has a separate method for parsing.

        parseMicroPhase(doc, cur, numEntry, kineticData);
      }
      cur = cur->next;
    }

    /// Push a copy of the isKinetic vector to the ChemicalSystem

    chemSys_->setIsKinetic(isKinetic_);

    xmlFreeDoc(doc);
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

void KineticController::parseMicroPhase(xmlDocPtr doc, xmlNodePtr cur,
                                        int &numEntry,
                                        struct KineticData &kineticData) {
  xmlChar *key;
  int proposedgemphaseid, proposedDCid;
  int testgemid, testdcid;
  string testname;
  bool kineticfound = false;
  bool ispozz = false;
  bool isParrotKilloh = false;
  bool istherm = false;
  bool issol = false;

  initKineticData(kineticData);

  cur = cur->xmlChildrenNode;

  isKinetic_.push_back(false);

  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"thamesname"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string testname((char *)key);
      kineticData.name = testname;
      kineticData.microPhaseId = chemSys_->getMicroPhaseId(testname);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"kinetic_data"))) {
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

      parseKineticData(doc, cur, kineticData);
    }

    cur = cur->next;
  }

  if (kineticfound) {
    kineticData.scaledMass =
        chemSys_->getMicroPhaseMass(kineticData.microPhaseId);
    kineticData.temperature = temperature_;
    kineticData.reftemperature = refT_;
    makeModel(doc, cur, kineticData);
  }

  /// Some items should be added to vectors whether kinetically controlled or
  /// not

  name_.push_back(kineticData.name);
  microPhaseId_.push_back(kineticData.microPhaseId);
  initScaledMass_.push_back(0.0);
  scaledMass_.push_back(0.0);

  return;
}

void KineticController::parseKineticData(xmlDocPtr doc, xmlNodePtr cur,
                                         struct KineticData &kineticData) {
  bool typefound = false;
  xmlChar *key;
  cur = cur->xmlChildrenNode;

  try {
    while (cur != NULL) {
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"type"))) {
        key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        string st((char *)key);
        kineticData.type = st;
        if (kineticData.type == ParrotKillohType) {
          typefound = true;
          parseKineticDataForParrotKilloh(doc, cur, kineticData);
        } else if (kineticData.type == StandardType) {
          typefound = true;
          parseKineticDataForStandard(doc, cur, kineticData);
        } else if (kineticData.type == PozzolanicType) {
          typefound = true;
          parseKineticDataForPozzolanic(doc, cur, kineticData);
        } else {
          xmlFree(key);
          throw HandleException("KineticController", "parseKineticData", "type",
                                "Model type not found");
        }
        xmlFree(key);
      }
      cur = cur->next;
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
    xmlDocPtr doc, xmlNodePtr cur, struct KineticData &kineticData) {
  xmlChar *key;
  cur = cur->next;

  if (verbose_) {
    cout << "--->Parsing PK data for " << kineticData.name << endl;
    cout.flush();
  }
  while (cur != NULL) {

    // Specific surface area (m2/kg)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"specificSurfaceArea"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.specificSurfaceArea, st);
      xmlFree(key);
    }
    // Reference specific surface area (m2/kg)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"refSpecificSurfaceArea"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.refSpecificSurfaceArea, st);
      xmlFree(key);
    }
    // Parrot-Killoh k1 parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"k1"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.k1, st);
      xmlFree(key);
    }
    // Parrot-Killoh k2 parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"k2"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.k2, st);
      xmlFree(key);
    }
    // Parrot-Killoh k3 parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"k3"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.k3, st);
      xmlFree(key);
    }
    // Parrot-Killoh n1 parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"n1"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.n1, st);
      xmlFree(key);
    }
    // Parrot-Killoh n3 parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"n3"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.n3, st);
      xmlFree(key);
    }
    // Parrot-Killoh HLK parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"HLK"))) {
        key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
        string st((char *)key);
        from_string(kineticData.HLK, st);
        xmlFree(key);
    }
    // Parrot-Killoh critical DOH parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"critdoh"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.critDOH, st);
      xmlFree(key);
    }
    // Activation energy
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"activationEnergy"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.activationEnergy, st);
      xmlFree(key);
    }
    cur = cur->next;
  }

  return;
}

void KineticController::parseKineticDataForStandard(
    xmlDocPtr doc, xmlNodePtr cur, struct KineticData &kineticData) {
  xmlChar *key;
  cur = cur->next;

  if (verbose_) {
    cout << "--->Parsing standard kinetic data for " << kineticData.name
         << endl;
    cout.flush();
  }
  while (cur != NULL) {

    // Specific surface area (m2/kg)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"specificSurfaceArea"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.specificSurfaceArea, st);
      xmlFree(key);
    }

    // Reference specific surface area (m2/kg)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"refSpecificSurfaceArea"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.refSpecificSurfaceArea, st);
      xmlFree(key);
    }

    // Dissolution rate constant (mol/m2/s)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dissolutionRateConst"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.dissolutionRateConst, st);
      xmlFree(key);
    }
    // Number of DC units produced in dissociation reaction
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dissolvedUnits"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.dissolvedUnits, st);
      xmlFree(key);
    }
    // Exponent on  the saturation index in the rate equation
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"siexp"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.siexp, st);
      xmlFree(key);
    }
    // Exponent on  the driving force term in the rate equation
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dfexp"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.dfexp, st);
      xmlFree(key);
    }
    // Loss on ignition of the material
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"loi"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.loi, st);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"activationEnergy"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.activationEnergy, st);
      xmlFree(key);
    }
    cur = cur->next;
  }

  return;
}

void KineticController::parseKineticDataForPozzolanic(
    xmlDocPtr doc, xmlNodePtr cur, struct KineticData &kineticData) {
  xmlChar *key;
  cur = cur->next;

  if (verbose_) {
    cout << "--->Parsing pozzolanic data for " << kineticData.name << endl;
    cout.flush();
  }
  while (cur != NULL) {

    // Specific surface area (m2/kg)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"specificSurfaceArea"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.specificSurfaceArea, st);
      xmlFree(key);
    }

    // Reference specific surface area (m2/kg)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"refSpecificSurfaceArea"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.refSpecificSurfaceArea, st);
      xmlFree(key);
    }

    // Dissolution rate constant (mol/m2/s)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dissolutionRateConst"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.dissolutionRateConst, st);
      xmlFree(key);
    }
    // Early-age diffusion rate constant (mol/m2/s)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"diffusionRateConstEarly"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.diffusionRateConstEarly, st);
      xmlFree(key);
    }
    // Number of DC units produced in dissociation reaction
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dissolvedUnits"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.dissolvedUnits, st);
      xmlFree(key);
    }
    // Later-age diffusion rate constant (mol/m2/s)
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"diffusionRateConstLate"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.diffusionRateConstLate, st);
      xmlFree(key);
    }
    // Exponent on  the saturation index in the rate equation
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"siexp"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.siexp, st);
      xmlFree(key);
    }
    // Exponent on  the driving force term in the rate equation
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dfexp"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.dfexp, st);
      xmlFree(key);
    }
    // Exponent on  the degree of reaction term in the diffusion rate equation
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dorexp"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.dorexp, st);
      xmlFree(key);
    }
    // Exponent on  the hydroxy ion activity in the rate equation
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"ohexp"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.ohexp, st);
      xmlFree(key);
    }
    // SiO2 mass fraction in the material
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"sio2"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.sio2, st);
      xmlFree(key);
    }
    // Al2O3 mass fraction in the material
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"al2o3"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.al2o3, st);
      xmlFree(key);
    }
    // CaO mass fraction in the material
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"cao"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.cao, st);
      xmlFree(key);
    }
    // Loss on ignition of the material
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"loi"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.loi, st);
      xmlFree(key);
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"activationEnergy"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.activationEnergy, st);
      xmlFree(key);
    }
    cur = cur->next;
  }

  return;
}

void KineticController::calcPhaseMasses(void) {
  int microPhaseId;
  double pscaledMass = 0.0;

  for (int i = 0; i < microPhaseId_.size(); i++) {
    microPhaseId = microPhaseId_[i];
    if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
      pscaledMass = chemSys_->getMicroPhaseMass(microPhaseId);
      scaledMass_[i] = pscaledMass;
      initScaledMass_[i] = pscaledMass;

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

  for (int i = 0; i < microPhaseId_.size(); i++) {
    microPhaseId = microPhaseId_[i];
    if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
      totmass += chemSys_->getMicroPhaseMass(microPhaseId);
    }
  }

  return (totmass);
}

void KineticController::makeModel(xmlDocPtr doc, xmlNodePtr cur,
                                  struct KineticData &kineticData) {
  KineticModel *km = NULL;

  if (kineticData.type == ParrotKillohType) {
    // Read remaining Parrot and Killoh model parameters
    km = new ParrotKillohModel(chemSys_, lattice_, kineticData,
                               verbose_, warning_);
  } else if (kineticData.type == StandardType) {
    // Read remaining pozzolanic model parameters
    km = new StandardKineticModel(chemSys_, lattice_, kineticData,
                                  verbose_, warning_);
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

  for (int midx = 0; midx < phaseKineticModel_.size(); ++midx) {
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
      cout << "KineticController::setPozzEffectOnPK for model " << midx << " ("
           << phaseKineticModel_[midx]->getType() << ")" << endl;
      cout << "  Ref LOI = " << loi << endl;
      cout << "  LOI = " << refloi << endl;
      cout << "  Max LOI = " << refloi << endl;
      cout << "  SiO2 = " << sio2val << endl;
      cout << "  Ref SiO2 = " << refsio2val << endl;
      cout << "  BET = " << betval << endl;
      cout << "  Ref BET = " << refbetval << endl;
      cout << "  Pozz Effect = " << pozzeffect << endl;
      cout << "  Min Pozz Effect = " << minpozzeffect << endl;
      cout.flush();
    }
  }

  minpozzeffect *= (refloi / maxloi);

  /// The way this is set up, 0.0 <= refloi / maxloi <= 1.0
  for (int midx = 0; midx < phaseKineticModel_.size(); ++midx) {
    if (phaseKineticModel_[midx]->getType() == ParrotKillohType) {
      phaseKineticModel_[midx]->setPfk(minpozzeffect);
    }
  }

  return;
}

void KineticController::calculateKineticStep (const double timestep,
                                             const double temperature,
                                             int cyc) {
    ///
    /// Initialize local variables
    ///
    ///

    double T = temperature;
    double arrhenius = 1.0;
    int i;

    //double massDissolved = 0.0;
    cout << scientific << setprecision(15);
    ///
    /// Determine if this is a normal step or a necessary
    /// tweak from a failed GEM_run call
    ///

    bool doTweak = (chemSys_->getTimesGEMFailed() > 0) ? true : false;

    static double hyd_time = 0.0;
    if (!doTweak)
        hyd_time = hyd_time + timestep;

    if (verbose_) {
        cout << "KineticController::calculateKineticStep Hydration Time = "
             << hyd_time << endl;
        cout.flush();
    }

    vector<double> impurityRelease;
    impurityRelease.clear();
    impurityRelease.resize(chemSys_->getNumMicroImpurities(), 0.0);
    vector <int> impurityDCID;
    impurityDCID.clear();
    impurityDCID.push_back(chemSys_->getDCId("K2O"));
    impurityDCID.push_back(chemSys_->getDCId("Na2O"));
    impurityDCID.push_back(chemSys_->getDCId("Per"));//170
    impurityDCID.push_back(chemSys_->getDCId("SO3"));

    //cout << endl << "impurityDCID : " << endl;
    //for(i = 0; i < chemSys_->getNumMicroImpurities(); i++){
    //    cout << i << "\t" << impurityDCID[i] << endl; cout.flush();
    //}
    //cout << endl ;

    //double molarMass;
    //vector<double> ICMoles, solutICMoles, DCMoles, GEMPhaseMoles;
    double totMassImpurity, massImpurity;

    int DCId;
    int pKMsize = phaseKineticModel_.size();
    static vector <double> scaledMass, massDissolved, keepNumDCMoles;
    static vector <int> phaseDissolvedId;
    double numDCMolesDissolved;

    try {

        //DCMoles_local.resize(DCNum_,0.0);


        cout << endl << "KineticController::calculateKineticStep     hyd_time = " << hyd_time << "\tcyc = " << cyc << endl;
        //cout << endl << "hydT     IC/DC/GEM/SOL - ini" << endl;

        for (int i = 0; i < DCNum_; i++) {
            DCMoles_[i] = chemSys_->getDCMoles(i);
            //cout <<"hydT   " << i << "\tDCName: " << chemSys_->getDCName(i)
            //     << "\tDCId: " << chemSys_->getDCId(chemSys_->getDCName(i)) << "\t" << DCMoles_[i] << endl; cout.flush();
        }

        for (int i = 0; i < ICNum_; i++) {
            ICMoles_[i] = 0.0;
        }

        if (hyd_time < leachTime_ && hyd_time < sulfateAttackTime_) {

            if (!doTweak) {
                // @todo BULLARD PLACEHOLDER
                // Still need to implement constant gas phase composition
                // Will involve equilibrating gas with aqueous solution
                //
                // First step each iteration is to equilibrate gas phase
                // with the electrolyte, while forbidding anything new
                // from precipitating.

                //vector<double> impurityRelease;= getDegreeOfReaction();
                //impurityRelease.clear();
                //impurityRelease.resize(chemSys_->getNumMicroImpurities(), 0.0);

                // RH factor is the same for all clinker phases
                //double vfvoid = lattice_->getVolumefraction(VOIDID);
                //double vfh2o = lattice_->getVolumefraction(ELECTROLYTEID);

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

                double critporediam = lattice_->getLargestSaturatedPore(); // in nm
                critporediam *= 1.0e-9;                                    // in m
                double rh = exp(-6.23527e-7 / critporediam / T);
                rh = rh > 0.55 ? rh : 0.551;

                /// Assume a zero contact angle for now.
                /// @todo revisit the contact angle issue

                /// Loop over all kinetic models
                //vector <double> rate;
                //rate.resize(pKMsize,0.0);
                scaledMass.resize(pKMsize,0.0);
                massDissolved.resize(pKMsize,0.0);
                phaseDissolvedId.resize(pKMsize,0.0);
                keepNumDCMoles.resize(pKMsize,0.0);
                double scMs, dsMs;

                //*******
                double totalDOR = (initScaledCementMass_  - chemSys_->getScaledCementMass()) / initScaledCementMass_;
                cout << "kineticController initScaledCementMass_/scaledCementMass : " << initScaledCementMass_
                     << " / " << chemSys_->getScaledCementMass() << endl;
                //*******

                for (int midx = 0; midx < pKMsize; ++midx) {
                    //phaseKineticModel_[midx]->calculateKineticStep(
                    //    timestep, temperature, isFirst, rh, dICMoles, dsolutICMoles,
                    //    tDCMoles, tGEMPhaseMoles);
                    DCId = phaseKineticModel_[midx]->getDCId();
                    phaseKineticModel_[midx]->calculateKineticStep(timestep, temperature, rh, scMs, dsMs, cyc, totalDOR);
                    scaledMass[midx] = scMs;
                    massDissolved[midx] = dsMs;
                    phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
                    chemSys_->setMicroPhaseMass(phaseDissolvedId[midx], scaledMass[midx]);
                    chemSys_->setMicroPhaseMassDissolved(phaseDissolvedId[midx], massDissolved[midx]);

                    if (verbose_) {
                        cout << "New scaled mass = "
                             << chemSys_->getMicroPhaseMass(phaseDissolvedId[midx])
                             << " and new volume = "
                             << chemSys_->getMicroPhaseVolume(phaseDissolvedId[midx]) << endl;
                        cout.flush();
                    }

                    //cout << "    midx/DCId/phaseDissolvedId : " << midx << " / " << DCId
                    //     << " / " << phaseDissolvedId[midx] << endl;
                    //cout << "      mass(%) K2O/Na2O/MgO/SO3 : " << chemSys_->getK2o(phaseDissolvedId[midx]) << " / "
                    //     << chemSys_->getNa2o(phaseDissolvedId[midx]) << //endl; cout.flush();
                    //     " / " << chemSys_->getMgo(phaseDissolvedId[midx]) <<
                    //     " / " << chemSys_->getSo3(phaseDissolvedId[midx]) << endl; cout.flush();

                    totMassImpurity = 0;

                    massImpurity = massDissolved[midx] * chemSys_->getK2o(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[0]] += massImpurity / chemSys_->getDCMolarMass("K2O");
                    //DCMoles_[impurityDCID[0]] += massDissolved[midx] * chemSys_->getK2o(phaseDissolvedId[midx]) / chemSys_->getDCMolarMass("K2O");

                    massImpurity = massDissolved[midx] * chemSys_->getNa2o(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[1]] += massImpurity / chemSys_->getDCMolarMass("Na2O");

                    massImpurity = massDissolved[midx] * chemSys_->getMgo(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[2]] += massImpurity / chemSys_->getDCMolarMass("Per"); //MgO

                    massImpurity = massDissolved[midx] * chemSys_->getSo3(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[3]] += massImpurity / chemSys_->getDCMolarMass("SO3");

                    //DCMoles_[DCId] += (massDissolved[midx] - totMassImpurity) / chemSys_->getDCMolarMass(DCId);
                    numDCMolesDissolved = (massDissolved[midx] - totMassImpurity) / chemSys_->getDCMolarMass(DCId);
                    keepNumDCMoles[midx] = DCMoles_[DCId] - numDCMolesDissolved;
                    chemSys_->setDCLowerLimit(DCId,keepNumDCMoles[midx]);

                    //check for tweak!!!

                }

               // cout << endl << " ******************** kinetic models ************************" << endl; //exit(0);
               // cout << " midx/name_/phaseDissolvedId[midx] :" << endl;
               // for (int midx = 0; midx < phaseKineticModel_.size(); ++midx) {
               //     cout << "  " << midx << "\tname: " << phaseKineticModel_[midx]->getName()
               //          << "\tmicroPhaseId_: " << phaseDissolvedId[midx]
               //          << "\tscldMass : " << scaledMass[midx] << "\tmassDiss : "
               //          << massDissolved[midx] << "\tDCId_loc : " << phaseKineticModel_[midx]->getDCId() << endl;
               // }
               //cout << endl << " ******************** kinetic models end ************************" << endl;exit(1);

            } else {
                //tweak!!!

                ///
                /// We will just tweak the icmoles a bit to try to
                /// cure a previous failed convergence with GEM_run
                ///

                //for (int ii = 0; ii < ICNum_; ii++) {
                //    if (ICName_[ii] != "H" && ICName_[ii] != "O" && ICName_[ii] != "Zz") {
                //        cout << "    " << ii << "    ICMoles_ini = " << ICMoles_[ii];
                //        ICMoles_[ii] *= 1.01;
                //        cout << "    ICMoles_fin = " << ICMoles[ii] << endl;
                //    }
                //}

                // correction on scaledMass/massdissolved of each kinetic controlled DC???

                double VMD = 0.01; // varMassDissolved
                cout << endl << "tweak for cyc = " << cyc << "   varMassDissolved = " << VMD << endl;
                cout << endl << "tweak before:" << endl;
                for (int midx = 0; midx < pKMsize; ++midx) {
                    cout << "   midx = " << midx << "     scaledMass[midx] = "
                         << scaledMass[midx] << "     massDissolved[midx] = "
                         << massDissolved[midx] << endl;
                }
                double deltaMassDissolved;

                for (int midx = 0; midx < pKMsize; ++midx) {
                    deltaMassDissolved = VMD * massDissolved[midx];
                    massDissolved[midx] += deltaMassDissolved;
                    scaledMass[midx] -= deltaMassDissolved;
                    phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
                    chemSys_->setMicroPhaseMass(phaseDissolvedId[midx], scaledMass[midx]);
                    chemSys_->setMicroPhaseMassDissolved(phaseDissolvedId[midx], massDissolved[midx]);

                    if (verbose_) {
                        cout << "New scaled mass = "
                             << chemSys_->getMicroPhaseMass(phaseDissolvedId[midx])
                             << " and new volume = "
                             << chemSys_->getMicroPhaseVolume(phaseDissolvedId[midx]) << endl;
                        cout.flush();
                    }
                    DCId = phaseKineticModel_[midx]->getDCId();
                    //cout << "    midx/pKMsize : " << midx << "/" << pKMsize
                    //     << " , " << chemSys_->getK2o(phaseDissolvedId[midx])
                    //     << " , " << chemSys_->getNa2o(phaseDissolvedId[midx])
                    //     << " , " << chemSys_->getMgo(phaseDissolvedId[midx])
                    //     << " , " << chemSys_->getSo3(phaseDissolvedId[midx]) << endl;

                    totMassImpurity = 0;

                    massImpurity = deltaMassDissolved * chemSys_->getK2o(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[0]] += massImpurity / chemSys_->getDCMolarMass("K2O");

                    massImpurity = deltaMassDissolved * chemSys_->getNa2o(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[1]] += massImpurity / chemSys_->getDCMolarMass("Na2O");

                    massImpurity = deltaMassDissolved * chemSys_->getMgo(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[2]] += massImpurity / chemSys_->getDCMolarMass("Per");//MgO

                    massImpurity = deltaMassDissolved * chemSys_->getSo3(phaseDissolvedId[midx]);
                    totMassImpurity += massImpurity;
                    DCMoles_[impurityDCID[3]] += massImpurity / chemSys_->getDCMolarMass("SO3");

                    //DCMoles_[DCId] += (deltaMassDissolved - totMassImpurity) / chemSys_->getDCMolarMass(DCId);
                    numDCMolesDissolved = (deltaMassDissolved - totMassImpurity) / chemSys_->getDCMolarMass(DCId);
                    keepNumDCMoles[midx] -= numDCMolesDissolved;
                    chemSys_->setDCLowerLimit(DCId,keepNumDCMoles[midx]);

                    //DCMoles_[DCId] += deltaMassDissolved / chemSys_->getDCMolarMass(DCId);
                    //DCMoles_[impurityDCID[0]] += deltaMassDissolved * chemSys_->getK2o(phaseDissolvedId[midx]) / chemSys_->getDCMolarMass("K2O");
                    //DCMoles_[impurityDCID[1]] += deltaMassDissolved * chemSys_->getNa2o(phaseDissolvedId[midx]) / chemSys_->getDCMolarMass("Na2O");
                    //DCMoles_[impurityDCID[2]] += deltaMassDissolved * chemSys_->getMgo(phaseDissolvedId[midx]) / chemSys_->getDCMolarMass("Per");//MgO
                    //DCMoles_[impurityDCID[3]] += deltaMassDissolved * chemSys_->getSo3(phaseDissolvedId[midx]) / chemSys_->getDCMolarMass("SO3");

                }

                cout << endl << "tweak after:" << endl;
                for (int midx = 0; midx < pKMsize; ++midx) {
                    cout << "   midx = " << midx << "     scaledMass[midx] = " << scaledMass[midx]
                         << "     massDissolved[midx] = " << massDissolved[midx] << endl;
                }

            }

            //cout << endl << "ICMoles_ in chemSys_ before transfer from kinetic in cyc = " << cyc << endl;
            //chemSys_->writeICMoles();

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
        EOBException ex("KineticController", "calculateKineticStep", oor.what(), 0, 0);
        ex.printException();
        exit(1);
    }

    for (i = 0; i < DCNum_; i++) {
        //cout << " " << i << "\t" << DCName_[i] << ": " << DCMoles_[i] << " mol" << endl;
        chemSys_->setDCMoles(i, DCMoles_[i]);
        //cout << "          " << DCName_[i] << ": " << chemSys_->getDCMoles(i) << " mol" << endl;
    }

    cout << "end calculateKineticStep - cyc = " << cyc << endl << endl;
    cout.flush();
    //exit(0);
    //if(cyc == 1){cout << "stop calculateKineticStep after cyc = " << cyc << endl;exit(0);}

    return;
}

