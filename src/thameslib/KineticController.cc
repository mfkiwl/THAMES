/**
@file  KineticController.cc
@brief Method definitions for the KineticController class.

*/
#include "KineticController.h"
#include "global.h"

KineticController::KineticController() {
  temperature_ = 293.15;

  // default temperature (K)
  refT_ = 293.15;

  ///
  /// Clear out the vectors so they can be populated with values from the
  ///

  numPhases_ = 0;
  chemSys_ = NULL;
  solut_ = NULL;
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

KineticController::KineticController(ChemicalSystem *cs, Solution *solut,
                                     Lattice *lattice, const string &fileName,
                                     const bool verbose, const bool warning)
    : chemSys_(cs), solut_(solut), lattice_(lattice) {
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
  ICName_ = chemSys_->getICName();
  DCName_ = chemSys_->getDCName();
  GEMPhaseNum_ = chemSys_->getNumGEMPhases();

  calcPhaseMasses();

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

  /// All kinetic components have been parsed now.
  /// Count the number of each type of kinetic model;

  int numpk = 0;
  int numpozz = 0;
  int numstd = 0;

  for (int i = 0; i < phaseKineticModel_.size(); ++i) {
    if (phaseKineticModel_[i]->getType() == ParrotKillohType) {
      numpk++;
    } else if (phaseKineticModel_[i]->getType() == PozzolanicType) {
      numpozz++;
    } else if (phaseKineticModel_[i]->getType() == StandardType) {
      numstd++;
    }
  }

  setNumPK(numpk);
  setNumPozzolanic(numpozz);
  setNumStandard(numstd);

  /// Next, this block tries
  /// to handle pozzolanic effects (loi, SiO2 content, etc.) on any other
  /// kinetic phases

  if (getNumPozzolanic() > 0) {
    setPozzEffectOnPK();
  }

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
    // Parrot-Killoh critical DOR parameter
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"critdoh"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      string st((char *)key);
      from_string(kineticData.critDOR, st);
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
    km = new ParrotKillohModel(chemSys_, solut_, lattice_, kineticData,
                               verbose_, warning_);
  } else if (kineticData.type == StandardType) {
    // Read remaining pozzolanic model parameters
    km = new StandardKineticModel(chemSys_, solut_, lattice_, kineticData,
                                  verbose_, warning_);
  } else if (kineticData.type == PozzolanicType) {
    // Read remaining pozzolanic model parameters
    km = new PozzolanicModel(chemSys_, solut_, lattice_, kineticData, verbose_,
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
  double minpozzeffect = 1000.0;
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

void KineticController::calculateDissolutionEvents(const double timestep,
                                                   const double temperature,
                                                   bool isFirst) {
  ///
  /// Initialize local variables
  ///

  double T = temperature;
  double arrhenius = 1.0;

  double rate = 1.0e-10; // Selected rate

  double massDissolved = 0.0;

  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  bool doTweak = (chemSys_->getTimesGEMFailed() > 0) ? true : false;

  static double hyd_time = 0.0;
  if (!doTweak)
    hyd_time = hyd_time + timestep;

  if (verbose_) {
    cout << "KineticController::calculateDissolutionEvents Hydration Time = "
         << hyd_time << endl;
    cout.flush();
  }

  try {
    static int conc_index = 0;
    int microPhaseId, DCId, ICId;
    double molarMass;
    vector<double> ICMoles, solutICMoles, DCMoles, GEMPhaseMoles;
    vector<double> dICMoles, dsolutICMoles, tDCMoles, tGEMPhaseMoles;
    ICMoles.clear();
    ICMoles.resize(ICNum_, 0.0);
    solutICMoles.clear();
    solutICMoles.resize(ICNum_, 0.0);
    DCMoles.clear();
    DCMoles.resize(DCNum_, 0.0);
    GEMPhaseMoles.clear();
    GEMPhaseMoles.resize(GEMPhaseNum_, 0.0);
    dICMoles.clear();
    dICMoles.resize(ICNum_, 0.0);
    dsolutICMoles.clear();
    dsolutICMoles.resize(ICNum_, 0.0);
    tDCMoles.clear();
    tDCMoles.resize(DCNum_, 0.0);
    tGEMPhaseMoles.clear();
    tGEMPhaseMoles.resize(GEMPhaseNum_, 0.0);
    string icn;

    // Populate IC moles before kinetic step
    for (int i = 0; i < ICNum_; i++) {
      ICMoles[i] = chemSys_->getICMoles(i);
      if (isFirst) {
        ICMoles[i] = 1.0e-9;
      } else {
        ICMoles[i] = chemSys_->getICMoles(i);
      }
    }

    // Populate DC moles before kinetic step
    for (int i = 0; i < DCNum_; i++) {
      DCMoles[i] = chemSys_->getDCMoles(i);
    }

    // Populate GEM phase moles before kinetic step
    for (int i = 0; i < GEMPhaseNum_; i++) {
      GEMPhaseMoles[i] = chemSys_->getGEMPhaseMoles(i);
    }

    solutICMoles = chemSys_->getSolution();

    if (isFirst) { // Beginning of special first-time setup tasks

      // Determine initial total solid mass
      double solidMass = 0.0;
      vector<double> initSolidMasses = getInitScaledMass();
      for (int i = 0; i < initSolidMasses.size(); ++i) {
        solidMass += initSolidMasses[i];
      }

      double waterMass = solidMass * lattice_->getWsratio();
      double waterMolarMass = chemSys_->getDCMolarMass(waterId_);
      double waterMoles = waterMass / waterMolarMass;
      microPhaseId = 0;
      double psMass, psVolume;
      double volume = 0.0;
      if (verbose_) {
        cout << "KineticController::calculateDissolutionEvents isFirst *** "
                "Initial "
                "solid mass = "
             << solidMass << endl;
        cout << "KineticController::calculateDissolutionEvents isFirst *** w/s "
                "ratio "
                "= "
             << lattice_->getWsratio() << endl;
        cout << "KineticController::calculateDissolutionEvents isFirst *** "
                "Initial "
                "water mass = "
             << waterMass << endl;
        cout << "KineticController::calculateDissolutionEvents isFirst *** "
                "Initial "
                "water moles = "
             << waterMoles << endl;
        cout.flush();

        for (int i = 0; i < microPhaseId_.size(); ++i) {
          cout << "KineticController::calculateDissolutionEvents "
               << "Initial MICROSTRUCTURE phase amount:" << endl;
          psMass = chemSys_->getMicroPhaseMass(microPhaseId_[i]);
          psVolume = chemSys_->getMicroPhaseVolume(microPhaseId_[i]);
          cout << "KineticController::calculateDissolutionEvents     "
               << chemSys_->getMicroPhaseName(microPhaseId_[i]) << " ("
               << microPhaseId_[i] << "): mass = " << psMass
               << ", vol = " << psVolume << endl;
          cout.flush();
        }
      }

      for (int i = 0; i < ICNum_; i++) {
        if (ICName_[i] == "H") {
          ICMoles[i] = (2.0 * waterMoles);
          solut_->setICMoles(i, ICMoles[i]);
          chemSys_->setDCMoles(waterId_, waterMoles);
        }
        if (ICName_[i] == "O") {
          ICMoles[i] = waterMoles;
          solut_->setICMoles(i, ICMoles[i]);
        }
      }

      double wmv = chemSys_->getNode()->DC_V0(chemSys_->getDCId(WaterDCName),
                                              chemSys_->getP(),
                                              chemSys_->getTemperature());
      chemSys_->setGEMPhaseMass(chemSys_->getGEMPhaseId(WaterGEMName),
                                waterMass);
      chemSys_->setGEMPhaseVolume(chemSys_->getGEMPhaseId(WaterGEMName),
                                  wmv / waterMoles);

      // Modify initial pore solution composition if desired
      // input units are mol/kgw
      // watermass is in units of grams, not kg

      double kgWaterMass = waterMass / 1000.0;

      map<int, double> isComp = chemSys_->getInitialSolutionComposition();
      map<int, double>::iterator p = isComp.begin();
      vector<double> ics;

      while (p != isComp.end()) {
        if (verbose_) {
          cout << "KineticController::calculateDissolutionEvents "
               << "modifying initial pore solution" << endl;
          cout.flush();
        }
        if (verbose_) {
          cout << "KineticController::calculateDissolutionEvents "
               << "--->Adding " << p->second << " mol/kgw of "
               << DCName_[p->first] << " to initial solution." << endl;
          cout.flush();
        }
        // Get the vector of IC compositions for this DC
        ics = chemSys_->getDCStoich(p->first);
        for (int ii = 0; ii < ics.size(); ++ii) {
          ICMoles[ii] += (p->second * ics[ii] * kgWaterMass);
        }
        p++;
      }

      // @todo BULLARD PLACEHOLDER
      // Still need to implement constant gas phase composition

    } // End of special first-time tasks

    if (hyd_time < leachTime_ && hyd_time < sulfateAttackTime_) {

      if (!doTweak) {

        // @todo BULLARD PLACEHOLDER
        // Still need to implement constant gas phase composition
        // Will involve equilibrating gas with aqueous solution
        //
        // First step each iteration is to equilibrate gas phase
        // with the electrolyte, while forbidding anything new
        // from precipitating.

        double vfvoid = lattice_->getVolumefraction(VOIDID);
        double vfh2o = lattice_->getVolumefraction(ELECTROLYTEID);

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

        /// Assume a zero contact angle for now.
        /// @todo revisit the contact angle issue

        /// Loop over all kinetic models

        double minmoles = 0.0;

        for (int midx = 0; midx < phaseKineticModel_.size(); ++midx) {

          // zero out the temporary chemical compositions before
          // passing them to the kinetic model

          fill(dICMoles.begin(), dICMoles.end(), minmoles);
          fill(dsolutICMoles.begin(), dsolutICMoles.end(), minmoles);
          fill(tDCMoles.begin(), tDCMoles.end(), minmoles);
          fill(tGEMPhaseMoles.begin(), tGEMPhaseMoles.end(), minmoles);

          phaseKineticModel_[midx]->calculateDissolutionEvent(
              timestep, temperature, isFirst, rh, dICMoles, dsolutICMoles,
              tDCMoles, tGEMPhaseMoles);

          for (int im = 0; im < ICMoles.size(); ++im) {

            ICMoles[im] += (dICMoles[im] - minmoles);
            if (verbose_) {
              cout << "Kinetically added " << dICMoles[im] << " moles of IC "
                   << chemSys_->getICName(im) << " to system" << endl;

              solutICMoles[im] += (dsolutICMoles[im] - minmoles);
            }

            for (int im = 0; im < DCMoles.size(); ++im) {
              DCMoles[im] += (tDCMoles[im] - minmoles);
            }

            for (int im = 0; im < GEMPhaseMoles.size(); ++im) {
              GEMPhaseMoles[im] += (tGEMPhaseMoles[im] - minmoles);
            }
          }
        }
      } else {

        ///
        /// We will just tweak the icmoles a bit to try to
        /// cure a previous failed convergence with GEM_run
        ///

        for (int ii = 0; ii < ICMoles.size(); ii++) {
          if (ICName_[ii] != "H" && ICName_[ii] != "O" && ICName_[ii] != "Zz") {
            ICMoles[ii] *= 1.01;
          }
        }
      }

      if (verbose_) {
        cout << "KineticController::calculateDissolutionEvents ICmoles after ";
        if (!doTweak) {
          cout << "dissolving:";
        } else {
          cout << "tweaking:";
        }
        cout << endl;
        for (int i = 0; i < ICNum_; i++) {
          cout << "    " << ICName_[i] << ": " << ICMoles[i] << " mol" << endl;
        }
        cout.flush();
      }

      if (doTweak) {
        for (int ii = 0; ii < ICMoles.size(); ii++) {
          chemSys_->setICMoles(ii, ICMoles[ii]);
        }
        return;
      }

      // Now correct for fixed solution composition, if any specified

      map<int, double> fsComp = chemSys_->getFixedSolutionComposition();
      if (fsComp.size() > 0) {

        // Get current mass of liquid water
        double kgWaterMass =
            chemSys_->getDCMoles(waterId_) * chemSys_->getDCMolarMass(waterId_);
        kgWaterMass *= 0.001; // molar mass is in g/mol

        map<int, double>::iterator p = fsComp.begin();
        vector<double> ics;

        while (p != fsComp.end()) {
          if (verbose_) {
            cout << "KineticController::calculateDissolutionEvents "
                 << "modifying pore solution" << endl;
            cout << "KineticController::calculateDissolutionEvents "
                 << "--->Adding " << p->second << " mol/kgw of "
                 << DCName_[p->first] << " to initial solution." << endl;
            cout.flush();
          }

          // Get current concentration of this DC
          double currentDCMoles = chemSys_->getDCMoles(p->first);
          // Convert it to molal concentration
          double currentDCConc = currentDCMoles / kgWaterMass;

          // Get difference in concentration from target value
          double diffConc = p->second - currentDCConc;
          // Convert the difference to mole difference
          double diffMoles = diffConc * kgWaterMass;

          // Now we know how many moles of this DC to add to the solution
          // Get the vector of IC compositions for this DC
          ics = chemSys_->getDCStoich(p->first);
          if (verbose_)
            cout << "Fixed BC effect:" << endl;
          for (int ii = 0; ii < ics.size(); ++ii) {
            ICMoles[ii] += (diffMoles * ics[ii]);
            if (verbose_) {
              cout << "    Add " << diffMoles * ics[ii] << " moles of IC "
                   << chemSys_->getICName(ii) << " due to BCs" << endl;
              cout << "    New IC concentration = " << ICMoles[ii] / kgWaterMass
                   << endl;
              cout.flush();
            }
          }
          p++;
        }
      }

      // Finally, set the new IC moles in the thermodynamic system
      for (int ii = 0; ii < ICMoles.size(); ii++) {
        chemSys_->setICMoles(ii, ICMoles[ii]);
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
    EOBException ex("KineticController", "calculateDissolutionEvents",
                    oor.what(), 0, 0);
    ex.printException();
    exit(1);
  }

  return;
}
