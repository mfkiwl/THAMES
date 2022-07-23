/**
@file  KineticModel.cc
@brief Method definitions for the KineticModel class.

*/
#include "KineticModel.h"

KineticModel::KineticModel ()
{
    ///
    /// Default value for w/c ratio in PK model is 0.45
    ///

    wcRatio_ = 0.45;

    ///
    /// Default value for Blaine fineness in PK model is 385 m<sup>2</sup>/kg
    ///

    blaine_ = 385.0;
    refBlaine_ = 385.0;     // reference Blaine fineness (m2/kg)

    ///
    /// Default temperature in the PK model is 20 C (or 293 K)
    ///

    temperature_ = 293.15;  // default temperature (K)
    refT_ = 293.15;         // default temperature (K)

    ///
    /// Clear out the vectors so they can be populated with values from the
    /// XML input file
    ///

    name_.clear();
    isKinetic_.clear();
    isThermo_.clear();
    isSoluble_.clear();
    microPhaseId_.clear();
    DCId_.clear();
    GEMPhaseId_.clear();
    RdICId_.clear();
    k1_.clear();
    k2_.clear();
    k3_.clear();
    n1_.clear();
    n3_.clear();
    Rd_.clear();
    activationEnergy_.clear();
    scaledMass_.clear();
    initScaledMass_.clear();
    critDOH_.clear();
    degreeOfHydration_.clear();
    NaTarget_.clear();
    KTarget_.clear();
    MgTarget_.clear();
    SO4Target_.clear();
    
    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10;

    return;
}

KineticModel::KineticModel (ChemicalSystem *cs,
                            Solution *solut,
                            Lattice *lattice,
                            const string &fileName,
                            const bool verbose,
                            const bool warning,
                            const bool debug)
:chemSys_(cs),solut_(solut),lattice_(lattice)
{
    const string PHASENUM = "phasenum";
    const string BEGINPHASE = "<phase>";
    const string ENDPHASE = "</phase>";
    const string PHASENAME = "name";
    const string TYPE = "type";
    const string GEMNAME = "gemname";
    const string DCNAME = "dcname";
    const string BLAINE = "Blaine";
    const string REFBLAINE = "refBlaine";
    const string WCRATIO = "wcRatio";
    const string TEMPERATURE = "Temperature";
    const string REFTEMPERATURE = "refTemperature";
    const string K1 = "K1";
    const string K2 = "K2";
    const string K3 = "K3";
    const string N1 = "N1";
    const string N3 = "N3";
    const string RD = "Rd";
    const string EA = "Ea";
    const string SCALEDMASS = "scaledMass";
    const string CRITDOH = "criticalDOH";
    
    // Set the verbose and warning flags
   
    verbose_ = verbose;
    warning_ = warning;
    debug_ = debug;

    ///
    /// Default value for Blaine fineness in PK model is 385 m<sup>2</sup>/kg
    ///

    blaine_ = 385.0;
    refBlaine_ = 385.0;
    blaineFactor_ = 1.0;    // ratio of blaine_ to refBlaine_

    ///
    /// Default value for w/c ratio in PK model is 0.45
    ///

    wcRatio_ = 0.45;

    ///
    /// Default initial solid mass is 100 g
    ///
    
    initSolidMass_ = 100.0;

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
    isKinetic_.clear();
    isThermo_.clear();
    isSoluble_.clear();
    microPhaseId_.clear();
    DCId_.clear();
    GEMPhaseId_.clear();
    RdICId_.clear();
    k1_.clear();
    k2_.clear();
    k3_.clear();
    n1_.clear();
    n3_.clear();
    Rd_.clear();
    activationEnergy_.clear();
    scaledMass_.clear();
    initScaledMass_.clear();
    critDOH_.clear();
    degreeOfHydration_.clear();
    NaTarget_.clear();
    KTarget_.clear();
    MgTarget_.clear();
    SO4Target_.clear();
    
    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
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
          if (verbose_) cout << "KineticModel data file is an XML file" <<endl;
          parseDoc(fileName);
      } else {
          throw FileException("KineticModel","KineticModel",fileName,
                            "NOT in XML format");
      }
    }
    catch (FileException fex) {
      fex.printException();
      exit(1);
    }

    int microPhaseId;
    if (verbose_) {
        cout << "RR Finished reading chemistry.xml file" << endl;
        for (int i = 0; i < microPhaseId_.size(); ++i) {
            microPhaseId = microPhaseId_[i];
            if (isKinetic(i)) {
                cout << "RR kinetic phase " << microPhaseId << endl;
                cout << "RR     name = " << chemSys_->getMicroPhaseName(microPhaseId)
                     << endl;
                cout << "RR     k1 =  " << k1_[i] << endl;
                cout << "RR     k2 =  " << k2_[i] << endl;
                cout << "RR     k3 =  " << k3_[i] << endl;
                cout << "RR     n1 =  " << n1_[i] << endl;
                cout << "RR     n3 =  " << n3_[i] << endl;
                cout << "RR     Ea =  " << activationEnergy_[i] << endl;
            }
        }
    }

    // JWB: We now know the identity of each phase in the microstructure and can
    // link it to a GEMS phase.  Now is the time to set the w/s ratio and
    // the initial mass fractions of the kinetic phases, rather than letting
    // the kinetic input file tell us.  This forces a direct link to the microstructure
    
    setInitialPhaseVolumeFractions();

    // JWB: Normalize mass fractions to total system mass and multiply by 100
    
    // Finally, set the initial total volume of the microstructure here
    
    double totmicvol = 0.0;
    for (int i = 0; i < microPhaseId_.size(); i++) {
        microPhaseId = microPhaseId_[i];
        if (microPhaseId != VOIDID) {
            totmicvol += chemSys_->getMicroPhaseVolume(microPhaseId);
        }
    }
    chemSys_->setMicroInitVolume(totmicvol);

    return;
}

void KineticModel::parseDoc (const string &docName)
{
    int numEntry = -1;   // Tracks number of solid phases
    int testgemid;

    ///
    /// The kineticData structure is used to temporarily hold parsed data
    /// for a given phase before the data are loaded permanently into class members.
    ///

    KineticData kineticData;

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
    ///

    try {
        string rxcsd = xsd_files_path;
        rxcsd+="/chemistry.xsd";
        if (verbose_) cout << "Chemistry xsd file is at " << rxcsd << endl;
        if(!is_xml_valid(doc,rxcsd.c_str())) {
            throw FileException("KineticModel","KineticModel",docName,
                                "xml NOT VALID");
        }

        if (doc == NULL ) {
            throw FileException("KineticModel","KineticModel",docName,
                            "xml NOT parsed successfully");
        }

        cur = xmlDocGetRootElement(doc);

        if (cur == NULL) {
            xmlFreeDoc(doc);
            throw FileException("KineticModel","KineticModel",docName,
                            "xml document is empty");
        }

        cur = cur->xmlChildrenNode;
        while (cur != NULL) {
            if ((!xmlStrcmp(cur->name, (const xmlChar *)"blaine"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(blaine_,st);
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"refblaine"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(refBlaine_,st);
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"wcRatio"))) {

                // Left in for backward compatibility, but this is
                // no longer used... w/c is calculated based on microstructure
                
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(wcRatio_,st);
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"temperature"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(temperature_,st);
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"reftemperature"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(refT_,st);
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"phase"))) {

                /// Each phase is a more complicated grouping of data that
                /// has a separate method for parsing.
                
                parsePhase(doc, cur, numEntry, kineticData);
            }
            cur = cur->next;
        }

        /// Push a copy of the isKinetic vector to the ChemicalSystem
        
        chemSys_->setIsKinetic(isKinetic_);

        xmlFreeDoc(doc);
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }
    return;
}

void KineticModel::parsePhase (xmlDocPtr doc,
                               xmlNodePtr cur,
                               int &numEntry,
                               KineticData &kineticData)
{
    xmlChar *key;
    int proposedgemphaseid,proposedDCid;
    int testgemid,testdcid;
    string testname;
    bool kineticfound = false;
    bool iskin = false;
    bool istherm = false;
    bool issol = false;

    initKineticData(kineticData);

    cur = cur->xmlChildrenNode;

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
            kineticData.GEMPhaseId =
                chemSys_->getMicroPhaseToGEMPhase(kineticData.microPhaseId,0);
            kineticData.DCId = chemSys_->getMicroPhaseToDC(kineticData.microPhaseId,0);

            ///
            /// Kinetic data are grouped together,
            /// so there is a method written just for parsing that grouping
            ///

            parseKineticData(doc, cur, kineticData);
        }
        cur = cur->next;
    }

    if (kineticfound) {

        if (kineticData.type == "kinetic") {
            if (verbose_) {
                cout << "QQ Kinetic Phase " << kineticData.name << ", id = "
                     << kineticData.microPhaseId << endl;
                cout << "QQ   (k1,k2,k3) =  " << kineticData.k1 << "," << kineticData.k2
                     << "," << kineticData.k3 << endl;
                cout << "QQ   (n1,n3) =  " << kineticData.n1 << "," << kineticData.n3 << endl;
                cout << "QQ   Ea = " << kineticData.Ea << endl;
            }
            iskin = true;
            istherm = false;
            issol = false;

        } else if (kineticData.type == "soluble") {

            iskin = false;
            istherm = true;
            issol = true;

        } else {

            iskin = false;
            istherm = true;
            issol = false;
        }
    } else {
        iskin = istherm = issol = false;
    }
    
    isKinetic_.push_back(iskin);
    isThermo_.push_back(istherm);
    isSoluble_.push_back(issol);
    initScaledMass_.push_back(0.0);
    scaledMass_.push_back(0.0);
    name_.push_back(kineticData.name);
    microPhaseId_.push_back(kineticData.microPhaseId);
    DCId_.push_back(kineticData.DCId);
    GEMPhaseId_.push_back(kineticData.GEMPhaseId);
    RdICId_.push_back(kineticData.RdId);
    Rd_.push_back(kineticData.RdVal);

    /// @note k1, k2, k3, n1, n3, and critDOH are all
    /// specific to PK model.
    
    /// Not knowing the order in which the clinker
    /// phases will appear, we must figure that out
    /// now.
    
    k1_.push_back(kineticData.k1);
    k2_.push_back(kineticData.k2);
    k3_.push_back(kineticData.k3);
    n1_.push_back(kineticData.n1);
    n3_.push_back(kineticData.n3);
    activationEnergy_.push_back(kineticData.Ea);
    critDOH_.push_back(kineticData.critDOH);

    return;
}

void KineticModel::parseKineticData (xmlDocPtr doc,
                                     xmlNodePtr cur,
                                     KineticData &kineticData)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    /// @note Everything in here pretty much depends on the PK
    /// model for PC clinker phase dissolution.
    ///
    /// @todo Replace PK model with more generic kinetic model
    /// for dissolution of each cement clinker phase and any
    /// others that may be kinetically controlled.
    
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"type"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            kineticData.type = st;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"scaledMass"))) {

            // Left in for backward compatibiltiy, but scaledMass
            // is no longer used here... it is calculated directly
            // from the microstructure
            
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.scaledMass,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k1"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.k1,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k2"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.k2,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k3"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.k3,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"n1"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.n1,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"n3"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.n3,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Ea"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.Ea,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"critdoh"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.critDOH,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Rd"))) {

            ///
            /// The data about partitioning of impurities among the clinker
            /// phases are grouped within a complex field in the input XML
            /// file, so we have a special method to parse it.
            ///

            parseRdData(doc, cur, kineticData);
        }
        cur = cur->next;
    }

    return;
}

void KineticModel::parseRdData(xmlDocPtr doc,
                               xmlNodePtr cur,
                               KineticData &kineticData) 
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;
    int RdId;
    double RdVal;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Rdelement"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            RdId = chemSys_->getICId(st);
            kineticData.RdId.push_back(RdId);
            xmlFree(key);
        }

        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Rdvalue"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(RdVal,st);
            kineticData.RdVal.push_back(RdVal);
            xmlFree(key);
        }
        cur = cur->next;
    }
}

void KineticModel::setInitialPhaseVolumeFractions()
{
    vector<double> microPhaseMass;
    microPhaseMass.clear();
    microPhaseMass.resize(microPhaseId_.size(),0.0);

    double volumeFraction,solidMass;
    double waterMass,molarMass,molarVolume,density;
    int DCId,microPhaseId;
    solidMass = 0.0;
    for (int i = 0; i < microPhaseId_.size(); ++i) {
        microPhaseId = microPhaseId_[i];
        if (microPhaseId != VOIDID) {
            if (verbose_) {
                cout << "  Found microstructure id " << microPhaseId;
                cout.flush();
                cout << " belongs to phase " << chemSys_->getMicroPhaseName(microPhaseId);
                cout.flush();
                volumeFraction = lattice_->getVolumefraction(microPhaseId);
                cout << "    Volume fraction = " << volumeFraction << endl;
                DCId = chemSys_->getMicroPhaseToDC(microPhaseId,0);
                cout << "    DC id number = " << DCId << endl;
                cout.flush();
                molarMass = chemSys_->getDCMolarMass(DCId);      // g/mol
                cout << "    Molar mass = " << molarMass << " g/mol" << endl;
                cout.flush();
                molarVolume = chemSys_->getDCMolarVolume(DCId);  // m3/mol
                cout << "    Molar volume = " << molarVolume << " m3/mol" << endl;
                cout.flush();
                density = molarMass / molarVolume / 1.0e6;          // g/cm3
                cout << "    Density = " << density << " g/cm3" << endl;
                cout.flush();
                microPhaseMass[microPhaseId] = volumeFraction * density;
                if (microPhaseId != WATERID) {
                    solidMass += microPhaseMass[microPhaseId];
                }
                cout << "    Mass = " << microPhaseMass[microPhaseId] << " g" << endl;
                cout.flush();
            } else {
                volumeFraction = lattice_->getVolumefraction(microPhaseId);
                DCId = chemSys_->getMicroPhaseToDC(microPhaseId,0);
                molarMass = chemSys_->getDCMolarMass(DCId);      // g/mol
                molarVolume = chemSys_->getDCMolarVolume(DCId);  // m3/mol
                density = molarMass / molarVolume / 1.0e6;          // g/cm3
                microPhaseMass[microPhaseId] = volumeFraction * density;
                if (microPhaseId != WATERID) {
                    solidMass += microPhaseMass[microPhaseId];
                }
            }
        }
    }

    // The water/solids mass ratio follows from that
   
    if (verbose_) {
        cout << "Input file thinks w/c = " << wcRatio_ << endl;
        cout.flush();
    }

    wcRatio_ = microPhaseMass[WATERID] / solidMass;

    if (verbose_) {
        cout << "Microstructure w/c = " << wcRatio_ << endl;
        cout.flush();
    }

    normalizePhaseMasses(microPhaseMass,solidMass);

    return;
}

void KineticModel::normalizePhaseMasses(vector<double> microPhaseMass,
                                        double solidMass)
{
    int microPhaseId;
    double pscaledMass = 0.0;

    for (int i = 0; i < microPhaseId_.size(); i++) {
        microPhaseId = microPhaseId_[i];
        if (microPhaseId == WATERID) {
            int waterId = chemSys_->getDCId("H2O@");
            double waterMolarMass = chemSys_->getDCMolarMass(waterId);
            pscaledMass = wcRatio_ * 100.0;  // Mass of solids scaled to 100 g now
            if (verbose_) {
                cout << "Setting DC moles of water to "
                     << (pscaledMass / waterMolarMass) << endl;
                cout.flush();
            }
            chemSys_->setDCMoles(waterId,(pscaledMass / waterMolarMass));
            if (verbose_) {
                cout << "Setting initial micphase mass and volume of "
                     << chemSys_->getMicroPhaseName(WATERID) << endl;
                cout.flush();
            }
            chemSys_->setMicroPhaseMass(WATERID,pscaledMass);
            chemSys_->setMicroPhaseMassDissolved(WATERID,0.0);

        } else if (microPhaseId != VOIDID) {

            pscaledMass = microPhaseMass[microPhaseId] * 100.0 / solidMass;
            if (isKinetic(i)) {
                scaledMass_[i] = pscaledMass;
                initScaledMass_[i] = pscaledMass;
                if (verbose_) {
                    cout << "Microstructure scaled mass of "
                         << chemSys_->getMicroPhaseName(microPhaseId)
                         << " (" << microPhaseId << ") = "
                         << microPhaseMass[microPhaseId] << " g out of "
                         << solidMass << " g total" << endl;
                    cout << "   Implies scaledMass_["
                         << i << "] = " << scaledMass_[i] << endl;
                }
            } else {

                if (isThermo(i)) {
                    scaledMass_[i] = pscaledMass;
                    initScaledMass_[i] = pscaledMass;
                    if (verbose_) {
                        cout << "Microstructure scaled mass of "
                             << chemSys_->getMicroPhaseName(microPhaseId)
                             << " (" << microPhaseId << ") = "
                             << microPhaseMass[microPhaseId] << " g out of "
                             << solidMass << " g total" << endl;
                        cout << "   Implies scaledMass_["
                             << i << "] = " << scaledMass_[i] << endl;
                        cout.flush();
                    }
                } else {
                    if (verbose_) {
                        cout << "Microstructure scaled mass of "
                             << chemSys_->getMicroPhaseName(microPhaseId)
                             << " (" << microPhaseId << ") = "
                             << microPhaseMass[microPhaseId] << " g out of "
                             << solidMass << " g total" << endl;
                        cout.flush();
                    }
                }
            }

            // Setting the phase mass will also automatically calculate the phase volume
            
            if (verbose_) {
                cout << "Setting initial micphase mass and volume of "
                     << chemSys_->getMicroPhaseName(microPhaseId) << endl;
                cout.flush();
            }
            chemSys_->setMicroPhaseMass(microPhaseId,pscaledMass);

            chemSys_->setMicroPhaseMassDissolved(microPhaseId,0.0);
        }
    }

    return;
}

void KineticModel::calculateKineticStep (const double timestep,
                                         const double temperature,
                                         bool isFirst)
{
    ///
    /// Initialize local variables
    ///

    double T = temperature;
    double wcFactor = 1.0;
    double DOH = 0.0;
    double arrhenius = 1.0;


    double ngrate = 1.0e-10;          // Nucleation and growth rate
    double hsrate = 1.0e-10;          // Hydration shell rate
    double diffrate = 1.0e-10;        // Diffusion rate
    double rate = 1.0e-10;            // Selected rate

    double newDOH = 0.0;              // Updated value of doh
    double massDissolved = 0.0;

    ///
    /// Determine if this is a normal step or a necessary
    /// tweak from a failed GEM_run call
    ///
    
    bool doTweak = (chemSys_->getTimesGEMFailed() > 0) ? true : false;

    static double hyd_time = 0.0;
    if (!doTweak) hyd_time = hyd_time + timestep;
    if (verbose_) cout << "hyd_time = " << hyd_time << endl;

    try {
        static int conc_index = 0;     
        int ICNum = chemSys_->getNumICs();
        int DCNum = chemSys_->getNumDCs();
        int numGEMPhases = chemSys_->getNumGEMPhases();
        int microPhaseId,DCId,ICId;
        double molarMass,Rd;
        vector<double> ICMoles,solutICMoles,DCMoles,GEMPhaseMoles;
        ICMoles.clear();
        ICMoles.resize(ICNum,0.0);
        solutICMoles.clear();
        solutICMoles.resize(ICNum,0.0);
        DCMoles.clear();
        DCMoles.resize(DCNum,0.0);
        GEMPhaseMoles.clear();
        GEMPhaseMoles.resize(numGEMPhases,0.0);
        string icn;

        vector<string> ICName;
        ICName.clear();
        ICName.resize(ICNum," ");

        int waterId = chemSys_->getDCId("H2O@");

        if (verbose_) {
            cout << "KineticModel::calculateKineticstep ICmoles before dissolving:" << endl;
            cout << "DC moles of water = " << chemSys_->getDCMoles(waterId);
            cout.flush();
            for (int i = 0; i < ICNum; i++) {
              if (isFirst) {
                  ICMoles[i] = 1.0e-9;
              } else {
                  ICMoles[i] = chemSys_->getICMoles(i);
              }
              ICName[i] = chemSys_->getICName(i);
              cout << "    " << ICName[i] << ": " << ICMoles[i] << " mol" << endl;
            }
        } else {
            for (int i = 0; i < ICNum; i++) {
              ICMoles[i] = chemSys_->getICMoles(i);
              if (isFirst) {
                  ICMoles[i] = 1.0e-9;
              } else {
                  ICMoles[i] = chemSys_->getICMoles(i);
              }
              ICName[i] = chemSys_->getICName(i);
            }
        }

        for (int i = 0; i < DCNum; i++) {
          DCMoles[i] = chemSys_->getDCMoles(i);
        }
        for (int i = 0; i < numGEMPhases; i++) {
          GEMPhaseMoles[i] = chemSys_->getGEMPhaseMoles(i);
        }

        solutICMoles = chemSys_->getSolution();

        if (isFirst) {  // Beginning of special first-time setup tasks
            
            // Set the proper amount of water for the total solid mass
            // and the water-solid ratio
            
            // Determine total intial solid mass
            double solidMass = 0.0;
            vector<double> initSolidMasses = getInitScaledMass();
            for (int i = 0; i < initSolidMasses.size(); i++) {
                solidMass += initSolidMasses[i];
            }
            double waterMass = solidMass * getWcRatio();
            double waterMolarMass = chemSys_->getDCMolarMass(waterId);
            double waterMoles = waterMass / waterMolarMass;
            if (verbose_) {
                cout << "*** Initial solid mass = " << solidMass << endl;
                cout << "*** w/s ratio = " << getWcRatio() << endl;
                cout << "*** Initial water mass = " << waterMass << endl;
                cout << "*** Initial water moles = " << waterMoles << endl;
                cout << "***" << endl;
            }

            if (waterMass <= 0.0) {
                throw FloatException("Controller","calculateState",
                                     "Divide by zero error");
            }

            // Print out the initial volumes of microstructure phases
           
            if (verbose_) {
                microPhaseId = 0;
                double psMass,psVolume;
                double volume = 0.0;
                cout << "Initial MICROSTRUCTURE phase amounts:" << endl;
                for (int i = 0; i < microPhaseId_.size(); ++i) {
                    microPhaseId = microPhaseId_[i];
                    if (microPhaseId != VOIDID) {
                        psVolume = chemSys_->getMicroPhaseVolume(microPhaseId);
                        volume += psVolume;
                    }
                }
                for (int i = 0; i < microPhaseId_.size(); ++i) {
                    microPhaseId = microPhaseId_[i];
                    if (microPhaseId != VOIDID) {
                        psMass = chemSys_->getMicroPhaseMass(microPhaseId);
                        psVolume = chemSys_->getMicroPhaseVolume(microPhaseId);
                        cout << "    " << chemSys_->getMicroPhaseName(microPhaseId)
                             << " (" << microPhaseId << "): mass = " << psMass
                             << ", vol = " << psVolume
                             << ", volfrac = " << (psVolume / volume) << endl;
                        cout.flush();
                    }
                }
            }

            for (int i = 0; i < ICMoles.size(); i++) {
                if (ICName[i] == "H") {
                    if (verbose_) {
                        cout << "previous IC moles for H is: " << ICMoles[i] << endl;
                    }
                    ICMoles[i] = (2.0 * waterMoles);
                    solut_->setICMoles(i,ICMoles[i]);
                    chemSys_->setDCMoles(waterId,waterMoles);
                    if (verbose_) {
                        cout << "new ICmoles for H is: " << ICMoles[i] << endl;
                    }
                }
                if (ICName[i] == "O") {
                    if (verbose_) cout << "previous IC moles for O is: " << ICMoles[i] << endl;
                    ICMoles[i] = waterMoles;
                    solut_->setICMoles(i,ICMoles[i]);
                    if (verbose_) cout << "new ICmoles for O is: " << ICMoles[i] << endl;
               }
            }
            double wmv = chemSys_->getNode()->DC_V0(chemSys_->getDCId("H2O@"),
                                                    chemSys_->getP(),
                                                    chemSys_->getTemperature());
            chemSys_->setGEMPhaseMass(chemSys_->getGEMPhaseId("aq_gen"),waterMass);
            chemSys_->setGEMPhaseVolume(chemSys_->getGEMPhaseId("aq_gen"),wmv/waterMoles);

            if (verbose_) {
                cout << "Looping over soluble minerals.  ";
                cout << "Here is the list of them:" << endl;
                for (int i = 0; i < microPhaseId_.size(); i++) {
                    if (isSoluble(i)) {
                        cout << name_[i] << endl;
                    }
                }
                cout.flush();
            }
            for (int i = 0; i < microPhaseId_.size(); i++) {
                if (isSoluble(i)) {
                    microPhaseId = microPhaseId_[i];
                    massDissolved = scaledMass_[i];
                    scaledMass_[i] = 0.0;
                    chemSys_->setMicroPhaseMass(microPhaseId,
                                          scaledMass_[microPhaseId]);
                    chemSys_->setMicroPhaseMassDissolved(microPhaseId,
                                       massDissolved);
                    if (verbose_) {
                        cout << name_[i] << ": Mass dissolved = "
                             << massDissolved << endl;
                        cout.flush();
                        cout << name_[i] << ": DC molar mass of "
                             << chemSys_->getDCName(DCId_[i])
                             << " = " << chemSys_->getDCMolarMass(DCId_[i])
                             << endl;
                        cout.flush();
                    }
                    for (int ii = 0; ii < ICMoles.size(); ii++) {
                        ICMoles[ii] += ((massDissolved
                               / chemSys_->getDCMolarMass(DCId_[i]))
                               * chemSys_->getDCStoich(DCId_[i],ii));
                        if (verbose_) {
                            cout << name_[i] << ":     Total dissolved "
                                 << ICName[ii] << " = " << ICMoles[ii] << " mol "
                                 << endl;
                        }
                    }
                    if (verbose_) cout.flush();

                } // END OF IF SOLUBLE
            }

            if (verbose_) {
                cout << "Looping over thermo minerals.  Here is the list of them:" << endl;
                for (int i = 0; i < microPhaseId_.size(); i++) {
                    if (isThermo(i)) {
                        cout << name_[i] << endl;
                    }
                }
                cout.flush();
            }
            for (int i = 0; i < microPhaseId_.size(); i++) {
                if (isThermo(i)) {
                    for (int ii = 0; ii < ICMoles.size(); ii++) {
                        ICMoles[ii] += ((scaledMass_[i]
                                    / chemSys_->getDCMolarMass(DCId_[i]))
                                    * chemSys_->getDCStoich(DCId_[i],ii));
                        if (verbose_) {
                            cout << name_[i] << ":     Total IC " << ICName[ii] << " = "
                                 << ICMoles[ii] << " mol " << endl;
                            cout.flush();
                        }
                    }
                
                } // END OF IF THERMO
            }

            // Modify initial pore solution composition if desired
            // input units are mol/kgw
            // watermass is in units of grams, not kg
            
            double kgWaterMass = waterMass / 1000.0;
            
            map<int,double> isComp = chemSys_->getInitialSolutionComposition();
            map<int,double>::iterator p = isComp.begin();

            while (p != isComp.end()) {
                if (verbose_) {
                    cout << "--->Adding " << p->second << " mol/kgw of "
                         << ICName[p->first] << " to initial solution." << endl;
                    cout.flush();
                }
                ICMoles[p->first] += (p->second * kgWaterMass);
                p++;
            }

        }   // End of special first-time tasks

        if (hyd_time < leachTime_ && hyd_time < sulfateAttackTime_) { 

          if (verbose_) {
             cout << "Looping over kinetically controlled phases.  " << endl;
             cout.flush();
          }

          vector<double> impurityRelease;
          impurityRelease.clear();
          impurityRelease.resize(chemSys_->getNumMicroImpurities(),0.0);

          if (verbose_) cout << "SS In KineticModel::calculateKineticStep" << endl;
          for (int i = 0; i < microPhaseId_.size(); i++) {
            microPhaseId = microPhaseId_[i];
            if (isKinetic(i)) {
                if (verbose_) {
                    cout << "SS kinetic phase " << microPhaseId << endl;
                    cout << "SS     name = " << chemSys_->getMicroPhaseName(microPhaseId)
                         << endl;
                    cout << "SS     k1 =  " << k1_[i] << endl;
                    cout << "SS     k2 =  " << k2_[i] << endl;
                    cout << "SS     k3 =  " << k3_[i] << endl;
                    cout << "SS     n1 =  " << n1_[i] << endl;
                    cout << "SS     n3 =  " << n3_[i] << endl;
                    cout << "SS     Ea =  " << activationEnergy_[i] << endl;
                }

                if (initScaledMass_[i] > 0.0) {
                  DOH = (initScaledMass_[i] - scaledMass_[i]) /
                        (initScaledMass_[i]);
                } else {
                  throw FloatException("KineticModel","calculateKineticStep",
                             "initScaledMass_ = 0.0");
                }

                wcFactor = 1.0 + (3.333 *
                       (pow(((critDOH_[i] * wcRatio_) - DOH),4.0)));

                arrhenius = exp((activationEnergy_[i]/GASCONSTANT)*((1.0/refT_) - (1.0/T)));
                // if (verbose_) {
                //     cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ" << endl;
                //     cout << "Calculating Arrhenius effect:" << endl;
                //     cout << "    GASCONSTANT = " << GASCONSTANT << endl;
                //     cout << "    activationEnergy_[" << i << "] = " << activationEnergy_[i] << endl;
                //     cout << "    refT_ = " << refT_ << endl;
                //     cout << "    T = " << T << endl;
                //     cout << "    arrhenius factor = " << arrhenius << endl;
                //     cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ" << endl;
                // }

                if (DOH < 1.0 && !doTweak) {
        
                    if (fabs(n1_[i]) > 0.0) {
                      ngrate = (k1_[i]/n1_[i]) * (1.0 - DOH)
                             * pow((-log(1.0 - DOH)),(1.0 - n1_[i]));
                      ngrate *= (blaineFactor_);  // only used for the N+G rate
              
                      if (ngrate < 1.0e-10) ngrate = 1.0e-10;
                    } else {
                      throw FloatException("KineticModel","calculateKineticStep",
                                     "n1_ = 0.0");
                    }
        
                    hsrate = k3_[i] * pow((1.0 - DOH),n3_[i]);
                    if (hsrate < 1.0e-10) hsrate = 1.0e-10;

                    diffrate = (k2_[i] * pow((1.0 - DOH),(2.0/3.0))) /
                               (1.0 - pow((1.0 - DOH),(1.0/3.0)));
                    if (diffrate < 1.0e-10) diffrate = 1.0e-10;


                    rate = (ngrate < hsrate) ? ngrate : hsrate;
                    if (diffrate < rate) rate = diffrate;
                    rate *= (wcFactor * arrhenius);
                    newDOH = DOH + (rate * timestep);
                    if (verbose_) {
                        cout << "PK model for " << getName(i)
                             << ", timestep " << timestep
                             << ", oldDOH = " << DOH << ", new DOH = "
                             << newDOH << endl;
                        cout.flush();
                    }

                    /// @note This where we can figure out the volume dissolved
                    /// and link it back to the current volume to see how many
                    /// voxels need to dissolve
                    ///

                    /// @note This all depends on concept of degree of
                    /// hydration as defined by the PK model 

                    /// @todo Make this independent of PK model
                    
                    scaledMass_[i] = initScaledMass_[i] * (1.0 - newDOH);
                    massDissolved = (newDOH - DOH) * initScaledMass_[i];
                    chemSys_->setMicroPhaseMass(microPhaseId,scaledMass_[i]);
                    chemSys_->setMicroPhaseMassDissolved(microPhaseId,massDissolved);
                    if (verbose_) {
                        cout << "Original scaled mass = " << initScaledMass_[i]
                             << " and new scaled mass = "
                             << chemSys_->getMicroPhaseMass(microPhaseId)
                             << " and new volume = "
                             << chemSys_->getMicroPhaseVolume(microPhaseId) << endl;
                    }

                    /// @note impurityRelease index values are assumed to
                    /// be uniquely associated with particular chemical
                    /// elements
                    
                    /// @todo Make this more general so that any indexing
                    /// can be used
                   
                    /// @todo Allow any IC element to be an impurity, not
                    /// necessarily just the ones hard-coded here
                   
                    impurityRelease[0] = (massDissolved *
                            chemSys_->getK2o(microPhaseId));
                    impurityRelease[1] = (massDissolved *
                            chemSys_->getNa2o(microPhaseId));
                    impurityRelease[2] = (massDissolved *
                            chemSys_->getMgo(microPhaseId));
                    impurityRelease[3] = (massDissolved *
                            chemSys_->getSo3(microPhaseId));

                    for (int ii = 0; ii < ICMoles.size(); ii++) {
                      ICMoles[ii] += ((massDissolved
                                  / chemSys_->getDCMolarMass(DCId_[i]))
                                  * chemSys_->getDCStoich(DCId_[i],ii));
                      if (ICName[ii] == "O") {
                        // Dissolved K2O in this phase
                        if (chemSys_->isIC("K")) {
                            icn = "K";
                            molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                            icn = "O";
                            molarMass += chemSys_->getICMolarMass(icn);
                            ICMoles[ii] += (impurityRelease[0]/molarMass);
                        }
                        // Dissolved Na2O in this phase
                        if (chemSys_->isIC("Na")) {
                            icn = "Na";
                            molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                            icn = "O";
                            molarMass += chemSys_->getICMolarMass(icn);
                            ICMoles[ii] += (impurityRelease[1]/molarMass);
                        }
                        // Dissolved MgO in this phase
                        if (chemSys_->isIC("Mg")) {
                            icn = "Mg";
                            molarMass = chemSys_->getICMolarMass(icn);
                            icn = "O";
                            molarMass += chemSys_->getICMolarMass(icn);
                            ICMoles[ii] += (impurityRelease[2]/molarMass);
                        }
                        // Dissolved SO3 in this phase
                         if (chemSys_->isIC("S")) {
                             icn = "S";
                             molarMass = chemSys_->getICMolarMass(icn);
                             icn = "O";
                             molarMass += (3.0 * chemSys_->getICMolarMass(icn));
                             ICMoles[ii] += (3.0 * (impurityRelease[3]/molarMass));
                         }
                       } else if (ICName[ii] == "S") {
                         // Dissolved SO3  in this phase
                         icn = "S";
                         molarMass = chemSys_->getICMolarMass(icn);
                         icn = "O";
                         molarMass += (3.0 * chemSys_->getICMolarMass(icn));
                         ICMoles[ii] += (impurityRelease[3]/molarMass);
                       } else if (ICName[ii] == "K") {
                         // Dissolved K2O in this phase
                         icn = "K";
                         molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                         icn = "O";
                         molarMass += chemSys_->getICMolarMass(icn);
                         ICMoles[ii] += (2.0 * (impurityRelease[0]/molarMass));
                       } else if (ICName[ii] == "Na") {
                         // Dissolved Na2O in this phase
                         icn = "Na";
                         molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                         icn = "O";
                         molarMass += chemSys_->getICMolarMass(icn);
                         ICMoles[ii] += (2.0 * (impurityRelease[1]/molarMass));
                       } else if (ICName[ii] == "Mg") {
                         // Dissolved MgO in this phase
                         icn = "Mg";
                         molarMass = chemSys_->getICMolarMass(icn);
                         icn = "O";
                         molarMass += chemSys_->getICMolarMass(icn);
                         ICMoles[ii] += (impurityRelease[2]/molarMass);
                       }
    
                     }

                } else if (DOH < 1.0) {

                    ///
                    /// We will just tweak the icmoles a bit to try to
                    /// cure a previous failed convergence with GEM_run
                    ///
                
                    for (int ii = 0; ii < ICMoles.size(); ii++) {
                        if (ICName[ii] != "H" && ICName[ii] != "O" && ICName[ii] != "Zz") {
                            ICMoles[ii] *= 1.01;
                        }
                    }

                } else {
                  throw DataException("KineticModel","calculateKineticStep",
                                    "DOH >= 1.0");
                }   
              }  // END OF IF KINETIC
          }
        }	

        if (verbose_) {
            if (!doTweak) {
                cout << "KineticModel::calculateKineticstep ICmoles after dissolving:" << endl;
            } else {
                cout << "KineticModel::calculateKineticstep ICmoles after tweaking:" << endl;
            }
            for (int i = 0; i < ICNum; i++) {
              cout << "    " << ICName[i] << ": " << ICMoles[i] << " mol" << endl;
            }
        }

        if (doTweak) {
            for (int ii = 0; ii < ICMoles.size(); ii++) {
                chemSys_->setICMoles(ii,ICMoles[ii]);
            }
            return;
        }

        vector<int> phaseList;
        double ICMol,ICStoich;
        double phMass,dphMass,partitionedMoles,denom;
        double partitionedMoles_K = 0.0, partitionedMoles_Na = 0.0;
        if (verbose_) {
            cout << "Aqueous mass = " << chemSys_->getGEMPhaseMass("aq_gen") << endl;
        }
        for (int i = 0; i < microPhaseId_.size(); i++) {
            microPhaseId = microPhaseId_[i];
            if (isThermo(i)) {
                if (verbose_) {
                    cout << "Partitioning alkali for phase " << name_[i] << ": ";
                }
                phaseList = chemSys_->getMicroPhaseToGEMPhase(microPhaseId);
                // List of all phases making up this microstructure phase
                // Compute the mass of each phase and sum to get phase mass
                phMass = dphMass = 0.0;
                for (int j = 0; j < phaseList.size(); j++) {
                    // Mass will be in g, not kg
                    if (verbose_) {
                        cout << endl;
                        cout << "phaseName: "
                             << chemSys_->getGEMPhaseName(phaseList[j]) << endl;
                        cout << "phaseMass: "
                             << chemSys_->getGEMPhaseMass(phaseList[j]) << endl;
                    }
                    dphMass += (double)(chemSys_->getGEMPhaseMass(phaseList[j])
                                      - chemSys_->getPrevGEMPhaseMass(phaseList[j]));
                    phMass += (double)(chemSys_->getGEMPhaseMass(phaseList[j]));
                }
                if (verbose_) {
                    cout << "Mass is " << phMass << "; mass increment is "
                         << dphMass << endl;
                    cout.flush();
                }
                for (int j = 0; j < RdICId_[i].size(); j++) {
                    ICId = RdICId_[i][j];
                    Rd = Rd_[i][j];
                    ICMol = chemSys_->getGEMPhaseStoich("aq_gen",ICId);
                    if (verbose_) {
                        cout << " Partitioning IC "
                             << chemSys_->getICName(ICId) << " with aqueous moles = ";
                        cout << ICMol << endl;
                        cout << " Partitioning IC " << chemSys_->getICName(ICId);
                        cout << " with actual aqueous moles = ";
                        cout << chemSys_->getGEMPhaseStoich(chemSys_->getGEMPhaseId("aq_gen"),ICId) << endl;
                        cout << "moles of 'O' in the solution: ";
                        cout << chemSys_->getGEMPhaseStoich(chemSys_->getGEMPhaseId("aq_gen"),
                                                         chemSys_->getICId("O")) << endl;
                        cout.flush();
                    }
                    if (dphMass > 0.0 && Rd > 0.0 && ICMol > 0.0) {
                        // phase volume will be in m3, so transform to cm3
                        denom = 1.0 + (1.0e6 * chemSys_->getGEMPhaseVolume("aq_gen")
                                               / dphMass / Rd);
                        partitionedMoles = (ICMol/denom);
            
                        // be careful whether partitionedmoles is the normalized ones?
                        ICMol -= (partitionedMoles);
                        if (ICId == chemSys_->getICId("Na")) {
                            partitionedMoles_Na = partitionedMoles * 
                                chemSys_->getGEMPhaseStoich(chemSys_->getGEMPhaseId("aq_gen"),
                                          chemSys_->getICId("O"));
                            if (verbose_) {
                                cout << "partitionedMoles for Na: "
                                     << partitionedMoles << endl;
                                cout << "actually removed Na: " << partitionedMoles_Na
                                     << " moles" << endl;
                            }
                            ICMoles[ICId] = ICMoles[ICId] - partitionedMoles_Na;
                        }
                        if (ICId == chemSys_->getICId("K")) {
                            partitionedMoles_K = partitionedMoles * 
                            chemSys_->getGEMPhaseStoich(chemSys_->getGEMPhaseId("aq_gen"),
                                             chemSys_->getICId("O"));
                            if (verbose_) {
                                cout << "partitionedMoles for K: "
                                     << partitionedMoles << endl;
                                cout << "actually removed K: " << partitionedMoles_K
                                     << " moles" << endl;
                            }
                            ICMoles[ICId] = ICMoles[ICId] - partitionedMoles_K;
                        }
                        if (ICMol < 1.0e-16) ICMol = 1.0e-16;
                        if (verbose_) {
                            cout << "; new conc = "
                                 << ICMol / chemSys_->getGEMPhaseMass("aq_gen")
                                  * 1000.0 << " molal" << endl;
                            cout.flush();
                        }
                    }
                }
            }  // END OF IF THERMO
        }  

        if (hyd_time >= leachTime_ && hyd_time < sulfateAttackTime_) {

            ///
            /// This is the block for simulating leaching kinetics
            ///

            if (verbose_) cout << " begin to leach at time = " << hyd_time << endl;
            double removed_K = 0.0, removed_Na = 0.0, removed_Ca = 0.0, removed_S = 0.0;
            double removed_O = 0.0, removed_H = 0.0;
            int KId = chemSys_->getDCId("K+");
            int NaId = chemSys_->getDCId("Na+");
            int CaId = chemSys_->getDCId("Ca+2");
            int SId = chemSys_->getDCId("SO4-2");
	
            // moles of K+, Na+, Ca+2 in solution before 
            double K1_aq, Na1_aq, Ca1_aq, S1_aq;
            // moles of K+, Na+ Ca+2 in solution after 
            double K2_aq, Na2_aq, Ca2_aq, S2_aq;
            double waterMass = 0.0;
            waterMass = chemSys_->getDCMoles(chemSys_->getDCId("H2O@"))
                            * chemSys_->getDCMolarMass("H2O@");
            //mass will be in g, not kg
            K1_aq = (chemSys_->getNode())->Get_cDC((long int) KId)
                     * waterMass * 1.0e-3; //mass will be in g, not kg
            Na1_aq = (chemSys_->getNode())->Get_cDC((long int) NaId)
                     * waterMass * 1.0e-3; //mass will be in g, not kg
            Ca1_aq = (chemSys_->getNode())->Get_cDC((long int) CaId)
                     * waterMass * 1.0e-3; //mass will be in g, not kg
            S1_aq = (chemSys_->getNode())->Get_cDC((long int) SId)
                     * waterMass * 1.0e-3;
           
            K2_aq = ((chemSys_->getNode())->Get_cDC((long int) KId) * 0.6 > 0.0005) ? 
        	         (chemSys_->getNode())->Get_cDC((long int) KId) * 0.6 * waterMass * 1.0e-3 : 0.0005 * waterMass * 1.0e-3;
            Na2_aq = ((chemSys_->getNode())->Get_cDC((long int) NaId) * 0.6 > 0.0005) ?
        	         (chemSys_->getNode())->Get_cDC((long int) NaId) * 0.6 * waterMass * 1.0e-3 : 0.0005 * waterMass * 1.0e-3;
            S2_aq = 0.0005 * waterMass * 1.0e-3;
            // Ca2_aq = 0.0005 * waterMass * 1.0e-3;

            double targetConc = 0.0005;		
            if (verbose_) cout << "set target calcium concentration to be: "
                               << targetConc << endl;
    
            if ((chemSys_->getNode())->Get_cDC((long int) CaId) > (targetConc + 0.0001)) {
                Ca2_aq = 1.0e-7 * waterMass * 1.0e-3;
            } else {
                Ca2_aq = targetConc * waterMass * 1.0e-3;
            }
	
            if (verbose_) cout << "molar mass of water is: "
                               << chemSys_->getDCMolarMass("H2O@") << endl;
	
            removed_K = K1_aq - K2_aq;
            removed_Na = Na1_aq - Na2_aq;
            removed_Ca = Ca1_aq - Ca2_aq;
            if (S1_aq > S2_aq){
                removed_S = S1_aq - S2_aq;
            } else {
                removed_S = 0.0;
            }

            for (int ii=0; ii < ICMoles.size(); ii++) {
                if (ICName[ii] == "K") {
                    ICMoles[ii] = ICMoles[ii] - removed_K;
                    if (verbose_) cout << "removed " << removed_K << "moles of K." << endl;
                }
                if (ICName[ii] == "Na") {
                    ICMoles[ii] = ICMoles[ii] - removed_Na;
             	    if (verbose_) cout << "removed " << removed_Na << "moles of Na." << endl;
                }
                if (ICName[ii] == "Ca") {
                    ICMoles[ii] = ICMoles[ii] - removed_Ca;
             	    if (verbose_) cout << "removed " << removed_Ca << "moles of Ca." << endl;
                }
                if (ICName[ii] == "S") {
                    ICMoles[ii] = ICMoles[ii] - removed_S;
                    if (verbose_) cout << "removed " << removed_S << "moles of S." << endl;
                }
                if (ICName[ii] == "O") {
                    removed_O = removed_Na + removed_K + 2 * removed_Ca + 4 * removed_S;
                    ICMoles[ii] = ICMoles[ii] - removed_O;
             	    if (verbose_) cout << "removed " << removed_O << "moles of O." << endl;
                }
                if (ICName[ii] == "H") {
                    removed_H = removed_Na + removed_K + 2 * removed_Ca;
                    ICMoles[ii] = ICMoles[ii] - removed_H;
             	    if (verbose_) cout << "removed " << removed_H << "moles of H." << endl;
                }
            }

        }   // End of block on leaching kinetics

        for (int ii = 0; ii < ICMoles.size(); ii++) {
            chemSys_->setICMoles(ii,ICMoles[ii]);
            if (verbose_) cout << "ICmoles of " << ICName[ii] << " is: " << ICMoles[ii] << endl;
        }

        if (hyd_time >= sulfateAttackTime_ && hyd_time < leachTime_) {
    
            ///
            /// This is the block for simulating sulfate attack kinetics
            ///

            double add_Na = 0.0;
            double add_S = 0.0;
            double add_O = 0.0;
            double add_H = 0.0;
            double add_OH = 0.0;
            double remove_K = 0.0;
            double remove_Ca = 0.0;
            double remove_O = 0.0;
            double remove_H = 0.0;
            double remove_OH = 0.0;
 
            int Sid = chemSys_->getDCId("SO4-2");
            int Naid = chemSys_->getDCId("Na+");
            int Kid = chemSys_->getDCId("K+");
            int Caid = chemSys_->getDCId("Ca+2");
            int Caid1 = chemSys_->getDCId("CaOH+");
            int Hid = chemSys_->getDCId("H+");
            int OHid = chemSys_->getDCId("OH-");


            // moles of components in solution before 
            double Na1_aq, K1_aq, Mg1_aq, S1_aq, Ca1_aq, H1_aq, OH1_aq;

            // moles of solution in solution after
            double Na2_aq, K2_aq, Mg2_aq, S2_aq, Ca2_aq, H2_aq, OH2_aq;

            double waterMass = 0.0;
            waterMass = chemSys_->getDCMoles(chemSys_->getDCId("H2O@"))
                       * chemSys_->getDCMolarMass("H2O@");


            double c1; //current concentration
            double c2; //new concentration

            int KId = chemSys_->getDCId("K+");
            int NaId = chemSys_->getDCId("Na+");
            int CaId = chemSys_->getDCId("Ca+2");
            int CaId1 = chemSys_->getDCId("CaOH+");
            int SId = chemSys_->getDCId("SO4-2");

            ifstream in("targetConc.dat");
            double concBuff = 0.0;
            in >> concBuff;
            double NaTarget = concBuff;
            if (verbose_) cout << "NaTarget = " << NaTarget << endl;
            in >> concBuff;
            double STarget = concBuff;
            if (verbose_) cout << "STarget = " << STarget << endl;
            in >> concBuff;
            double CaTarget = concBuff;
            if (verbose_) cout << "CaTarget = " << CaTarget << endl;
            in >> concBuff;
            double Ca1Target = concBuff;
            if (verbose_) cout << "Ca1Target = " << Ca1Target << endl;
            in >> concBuff;
            double KTarget = concBuff;
            if (verbose_) cout << "KTarget = " << KTarget << endl;
            in >> concBuff;

            in.close();    

            //double NaTarget = 0.15313262, STarget = 0.00196016;
            //double CaTarget = 0.00012733, Ca1Target = 0.00034906;
            //double KTarget = 0.4069377;
    
            for (int j = 0; j < 4; j++) { 
              switch(j) {
                case 0: // Na+
                  c1 = (chemSys_->getNode())->Get_cDC((long int) NaId);
                  c2 = NaTarget;
                  Na1_aq = c1 * waterMass * 1.0e-3;
                  Na2_aq = c2 * waterMass * 1.0e-3;
                  add_Na += (Na2_aq - Na1_aq);
                  break;

                case 1: // SO4-2
                  c1 = (chemSys_->getNode())->Get_cDC((long int) SId);
                  c2 = STarget;
                  S1_aq = c1 * waterMass * 1.0e-3;
                  S2_aq = c2 * waterMass * 1.0e-3;
                  add_S = S2_aq - S1_aq;
                  add_O += (add_S * 4.0);
                  break;

                case 2: // K+
                  c1 = (chemSys_->getNode())->Get_cDC((long int) KId);
                  c2 = KTarget;
                  if (c1 >= 0.0005) {
                    K1_aq = c1 * waterMass * 1.0e-3;
                    K2_aq = c2 * waterMass * 1.0e-3;
                    remove_K += (K1_aq - K2_aq);
                  } else {
                    remove_K = 0.0;
                  }
                  break;

                case 3: // Ca2+ and CaOH+  
                  c1 = ((chemSys_->getNode())->Get_cDC((long int) CaId)+ 
                           (chemSys_->getNode())->Get_cDC((long int) CaId1));
                  c2 = CaTarget + Ca1Target;
                  Ca1_aq = c1 * waterMass * 1.0e-3;
                  Ca2_aq = c2 * waterMass * 1.0e-3;
                  remove_Ca += (Ca1_aq - Ca2_aq);
                  remove_O += (remove_K + remove_Ca * 2);
                  remove_H += (remove_K + remove_Ca * 2);
                  break;

                /*
                case 4: // H+
                  c1 = ((chemSys_->getNode())->Get_cDC((long int)HId));        
                  c2 = H_target;
                  H1_aq = c1 * waterMass * 1.0e-3;
                  H2_aq = c2 * waterMass * 1.0e-3;
                  add_H += H2_aq - H1_aq;
                  break;
                */
                /*
                case 4: // OH-
                  c1 = ((chemSys_->getNode())->Get_cDC((long int)OHId));
                  c2 = OH_target;
                  OH1_aq = c1 * waterMass * 1.0e-3;
                  OH2_aq = c2 * waterMass * 1.0e-3;
                  add_OH = OH2_aq - OH1_aq;
                  add_O += add_OH;
                  add_H += add_OH;
                  cout << "added " << add_OH << " moles of OH-." << endl;
                  break;
                */
                default:
                  break;
              }
            } // end of test
    
            for (int ii=0; ii < ICMoles.size(); ii++) {
                if (ICName[ii] == "Na") {
                    ICMoles[ii] = ICMoles[ii] + add_Na;
                    solutICMoles[ii] = solutICMoles[ii] + add_Na;
                    if (verbose_) cout << "added " << add_Na << " moles of Na." << endl;
                }
                if (ICName[ii] == "S") {
                    ICMoles[ii] = ICMoles[ii] + add_S;
                    solutICMoles[ii] = solutICMoles[ii] + add_S;
                    if (verbose_) cout << "added " << add_S << " moles of S." << endl;
                }
                if (ICName[ii] == "K") {
                    ICMoles[ii] = ICMoles[ii] - remove_K;
                    solutICMoles[ii] = solutICMoles[ii] - remove_K;
                    if (verbose_) cout << "removed " << remove_K << " moles of K." << endl;
                }
                if (ICName[ii] == "Ca") {
                    ICMoles[ii] = ICMoles[ii] - remove_Ca;
                    solutICMoles[ii] = solutICMoles[ii] - remove_Ca;
                    if (verbose_) cout << "removed " << remove_Ca << " moles of Ca." << endl;
                }
                if (ICName[ii] == "O") {
                    ICMoles[ii] = ICMoles[ii] + add_O - remove_O;
                    solutICMoles[ii] = solutICMoles[ii] + add_O - remove_O;
                    if (verbose_) cout << "added " << (add_O - remove_O) << " moles of O." << endl;
                }
                if (ICName[ii] == "H") {
                    ICMoles[ii] = ICMoles[ii] + add_H - remove_H;
                    solutICMoles[ii] = solutICMoles[ii] + add_H - remove_H;
                    if (verbose_) cout << "added " << (add_H - remove_H) << " moles of H." << endl;
                }
            }
   
            for (int ii = 0; ii < ICMoles.size(); ii++) {
                chemSys_->setICMoles(ii,ICMoles[ii]);
                solut_->setICMoles(ii,solutICMoles[ii]);
            }

        }       // End of sulfate attack block
    }

    catch (EOBException eex) {
        eex.printException();
        exit(1);
    }
    catch (DataException dex) { 
        dex.printException();
        exit(1);
    }
    catch (FloatException fex) {
        fex.printException();
        exit(1);
    }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","calculateKineticStep",
                           oor.what(),0,0);
        ex.printException();
        exit(1);
    }

	
    return;
}
 
 
void KineticModel::calculatePhaseChange (const int microPhaseId,
                                         double k,
                                         double gamma,
                                         double timestep)
{
    ///
    /// Initialize local variables
    ///

    double DCChange = 0.0;
    vector<int> GEMPhaseIndex;
    GEMPhaseIndex.clear();
    GEMPhaseIndex = chemSys_->getMicroPhaseToGEMPhase(microPhaseId);

    ///
    /// Should the next block be only for verbose output
    ///

    if (verbose_) {
        cout << "Microphase " << chemSys_->getMicroPhaseName(microPhaseId)
             << " contains phases: "
             << endl;
        for (int i = 0; i < GEMPhaseIndex.size(); i++) {
            cout << chemSys_->getGEMPhaseName(i) << endl;
        }    
    }

    double GEMPhaseMoles = 0.0;
    double *GEMPhaseMolesArray = chemSys_->getGEMPhaseMoles();
    for (int i = 0; i < GEMPhaseIndex.size(); i++) {
      GEMPhaseMoles += GEMPhaseMolesArray[i];
    }

    vector<double> GEMPhaseFrac;
    GEMPhaseFrac.clear();
    GEMPhaseFrac.resize(GEMPhaseIndex.size(), 0.0);

    for (int i = 0; i < GEMPhaseIndex.size(); i++) {
      double A = 0.0; // A is the surface area
      GEMPhaseFrac[i] = chemSys_->getGEMPhaseMoles(GEMPhaseIndex[i]) / GEMPhaseMoles;
      A = lattice_->getSurfaceArea(microPhaseId) * GEMPhaseFrac[i];
      double SI = 0.0; // SI is the saturation index for the phase i
      SI = solut_->getSI(GEMPhaseIndex[i]);
      
      vector<int> DCIndex;
      DCIndex.clear();
      DCIndex = chemSys_->getGEMPhaseDCMembers(GEMPhaseIndex[i]);
      if (verbose_) cout << "Phase "
                         << chemSys_->getGEMPhaseName(GEMPhaseIndex[i]) 
                         << " contains DC members: " << endl;
      for (int j = 0; j < DCIndex.size(); j++) {
        chemSys_->getDCName(DCIndex[j]);
      }
      

      vector<double> DCFrac;
      DCFrac.clear();
      DCFrac.resize(DCIndex.size(), 0.0);
      double DCMoles = 0.0;
      for (int j = 0; j < DCIndex.size(); j++) {
        DCMoles += chemSys_->getDCMoles(DCIndex[j]);
      }
      
      double surfaceArea = 0.0;
      for (int j = 0; j < DCIndex.size(); j++) {
        DCFrac[j] = chemSys_->getDCMoles(DCIndex[j]) / DCMoles;
        surfaceArea = A * DCFrac[j];
        DCChange = k * surfaceArea * pow((SI - 1),gamma);

        ///
        /// Convert time step from days to seconds, because DCChange
        /// has units of mol/s
        ///

        DCChange = DCChange * timestep * 24.0 * 60.0 * 60.0; // DCChange unit: moles/s

        double newDCMoles = chemSys_->getDCMoles(DCIndex[j]) + DCChange;
        chemSys_->setDCUpperLimit(DCIndex[j],newDCMoles);
        chemSys_->setDCLowerLimit(DCIndex[j],newDCMoles);

      }
    }

    return;
}

void KineticModel::setKineticDCMoles ()
{
    if (verbose_) {
        cout << "    In setKineMicDCmoles..." << endl;
        cout.flush();
    }

    try {
        int waterId = chemSys_->getDCId("H2O@");
        double waterMoles = chemSys_->getDCMoles(waterId);
        double waterMolarMass = chemSys_->getDCMolarMass(waterId);
        double waterMass = waterMoles * waterMolarMass;
        if (verbose_) {
            cout << "        "
                 << chemSys_->getDCName(waterId) << ": Mass = " << waterMass
                 << ", Molar mass = " << waterMolarMass
                 << endl;
        }
        for (int i = 0; i < microPhaseId_.size(); i++) {
            if (isKinetic(i)) {
                if (chemSys_->getDCMolarMass(DCId_[i]) <= 0.0) {
                    throw FloatException("KineticModel","setKineticDCmoles",
                                         "Divide by zero error");
                }
                if (verbose_) {
                    cout << "        Clinker phase "
                         << name_[i] << ": Mass = " << scaledMass_[i]
                         << ", Molar mass = " << chemSys_->getDCMolarMass(DCId_[i])
                         << endl;
                }
                chemSys_->setDCMoles(DCId_[i],(scaledMass_[i]
                                 / chemSys_->getDCMolarMass(DCId_[i])));
            } // END OF IF KINETIC
        }
    }
    catch (EOBException eex) {
        eex.printException();
        exit(1);
    }
    catch (FloatException fex) {
        fex.printException();
        exit(1);
    }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","setKineticDCMoles",
                           oor.what(),0,0);
        ex.printException();
        exit(1);
    }
    return;
}

void KineticModel::zeroKineticDCMoles ()
{
    try {
        for (int i = 0; i < microPhaseId_.size(); i++) {
            if (isKinetic(i)) {
                chemSys_->setDCMoles(DCId_[i],0.0);
            }
        }
    }
    catch (EOBException eex) {
        eex.printException();
        exit(0);
    }
    return;
}
