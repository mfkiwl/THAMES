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
    isPK_.clear();
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

    pfk1_ = pfk2_ = pfk3_ = 1.0;  // Default PK pozzolanic factors for no pozzolans

    return;
}

KineticModel::KineticModel (ChemicalSystem *cs,
                            Solution *solut,
                            Lattice *lattice,
                            const string &fileName,
                            const bool verbose,
                            const bool warning)
:chemSys_(cs),solut_(solut),lattice_(lattice)
{
    // Set the verbose and warning flags
   
    #ifdef DEBUG
        verbose_ = true;
        warning_ = true;
        cout << "KineticModel::KineticModel Constructor" << endl;
        cout.flush();
    #else
        verbose_ = verbose;
        warning_ = warning;
    #endif

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
    isPK_.clear();
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
          #ifdef DEBUG
              cout << "KineticModel data file is an XML file" <<endl;
          #endif
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

    #ifdef DEBUG
        cout << "KineticModel::KineticModel Finished reading chemistry.xml file" << endl;
        for (int i = 0; i < microPhaseId_.size(); ++i) {
            microPhaseId = microPhaseId_[i];
            if (isKinetic(i)) {
                cout << "KineticModel::KineticModel kinetic phase " << microPhaseId << endl;
                cout << "KineticModel::KineticModel     name = " << chemSys_->getMicroPhaseName(microPhaseId)
                     << endl;
                cout << "KineticModel::KineticModel     k1 =  " << k1_[i] << endl;
                cout << "KineticModel::KineticModel     k2 =  " << k2_[i] << endl;
                cout << "KineticModel::KineticModel     k3 =  " << k3_[i] << endl;
                cout << "KineticModel::KineticModel     n1 =  " << n1_[i] << endl;
                cout << "KineticModel::KineticModel     n3 =  " << n3_[i] << endl;
                cout << "KineticModel::KineticModel     Ea =  " << activationEnergy_[i] << endl;
            }
        }
        cout.flush();
    #endif

    // The ChemicalSystem and Lattice objects were constructed before we get
    // here, so now we know nearly everything we need to know about the
    // microstructure.  The only thing left to do is to populate the scaled
    // phase masses and read the w/s ratio, which the kinetic model needs
    
    getPhaseMasses();

    pfk1_ = pfk2_ = pfk3_ = 1.0;   // Default PK factors for no pozzolanic material
                                   
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
    /// @note This block requires schema file to be local

    try {
        string rxcsd = "chemistry.xsd";
        #ifdef DEBUG
            cout << "KineticModel::parseDoc Chemistry xsd file is at "
                 << rxcsd << endl;
            cout.flush();
        #endif
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

        /// Set the Blaine correction factor

        if (refBlaine_ > 0.0) {
            blaineFactor_ = blaine_ / refBlaine_;
        } else {
            blaineFactor_ = 1.0;
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
    bool isPK = false;
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
            kineticData.DCId =
                chemSys_->getMicroPhaseDCMembers(kineticData.microPhaseId,0);

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
            #ifdef DEBUG
                cout << "KineticModel::parsePhase Kinetic Phase " << kineticData.name
                     << ", id = "
                     << kineticData.microPhaseId << endl;
                cout << "KineticModel::parsePhase   (k1,k2,k3) =  " << kineticData.k1
                     << "," << kineticData.k2
                     << "," << kineticData.k3 << endl;
                cout << "KineticModel::parsePhase   (n1,n3) =  " << kineticData.n1 << ","
                     << kineticData.n3 << endl;
                cout << "KineticModel::parsePhase   Ea = " << kineticData.Ea << endl;
                cout.flush();
            #endif
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

    // @todo BULLARD PLACEHOLDER
    // Special kluges here for silica fume, which we have to call Silica-amorph
    // for now to make it compatible with the user interface as of 2022 Dec 29
   
    string sfume("SFUME");
    string siamorph("Silica-amorph");
    if (kineticData.name == sfume || kineticData.name == siamorph) {
        kineticData.name = siamorph;
        isPK_.push_back(false);
    } else {
        isPK_.push_back(iskin);
    }
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

void KineticModel::getPhaseMasses(void)
{
    int microPhaseId;
    double pscaledMass = 0.0;

    wcRatio_ = lattice_->getWsratio();

    for (int i = 0; i < microPhaseId_.size(); i++) {
        microPhaseId = microPhaseId_[i];
        if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
            pscaledMass = chemSys_->getMicroPhaseMass(microPhaseId);
            scaledMass_[i] = pscaledMass;
            initScaledMass_[i] = pscaledMass;

            // Setting the phase mass will also automatically calculate the phase volume
            
            #ifdef DEBUG
                cout << "KineticModel::getPhaseMasses Kinetic model reads solid micphase mass of "
                     << chemSys_->getMicroPhaseName(microPhaseId)
                     << " as " << initScaledMass_[i] << endl;
                cout.flush();
            #endif
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
    double rhFactor = 1.0;
    double DOH = 0.0;
    double cDOH = 0.0;
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

    #ifdef DEBUG
        cout << "KineticModel::calculateKineticStep Hydration Time = "
             << hyd_time << endl;
        cout.flush();
    #endif

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

        #ifdef DEBUG
            cout << "KineticModel::calculateKineticStep ICmoles before dissolving:" << endl;
            cout << "KineticModel::calculateKineticStep DC moles of water = " << chemSys_->getDCMoles(waterId);
            cout.flush();
            for (int i = 0; i < ICNum; i++) {
              if (isFirst) {
                  ICMoles[i] = 1.0e-9;
              } else {
                  ICMoles[i] = chemSys_->getICMoles(i);
              }
              ICName[i] = chemSys_->getICName(i);
              cout << "KineticModel::calculateKineticStep     " << ICName[i] << ": " << ICMoles[i] << " mol" << endl;
            }
        #else 
            for (int i = 0; i < ICNum; i++) {
              ICMoles[i] = chemSys_->getICMoles(i);
              if (isFirst) {
                  ICMoles[i] = 1.0e-9;
              } else {
                  ICMoles[i] = chemSys_->getICMoles(i);
              }
              ICName[i] = chemSys_->getICName(i);
            }
        #endif

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
            
            // Determine total initial solid mass
            double solidMass = 0.0;
            vector<double> initSolidMasses = getInitScaledMass();
            for (int i = 0; i < initSolidMasses.size(); i++) {
                solidMass += initSolidMasses[i];
            }
            double waterMass = solidMass * getWcRatio();
            double waterMolarMass = chemSys_->getDCMolarMass(waterId);
            double waterMoles = waterMass / waterMolarMass;
            if (verbose_) {
                cout << "KineticModel::calculateKineticStep *** Initial solid mass = "
                     << solidMass << endl;
                cout << "KineticModel::calculateKineticStep *** w/s ratio = "
                     << getWcRatio() << endl;
                cout << "KineticModel::calculateKineticStep *** Initial water mass = "
                     << waterMass << endl;
                cout << "KineticModel::calculateKineticStep *** Initial water moles = "
                     << waterMoles << endl;
                cout << "KineticModel::calculateKineticStep ***" << endl;
                cout.flush();
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
                cout << "KineticModel::calculateKineticStep Initial MICROSTRUCTURE phase amounts:" << endl;
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
                        cout << "KineticModel::calculateKineticStep     " << chemSys_->getMicroPhaseName(microPhaseId)
                             << " (" << microPhaseId << "): mass = " << psMass
                             << ", vol = " << psVolume
                             << ", volfrac = " << (psVolume / volume) << endl;
                        cout.flush();
                    }
                }
            }

            for (int i = 0; i < ICMoles.size(); i++) {
                if (ICName[i] == "H") {
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Previous IC moles for H is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    #endif
                    ICMoles[i] = (2.0 * waterMoles);
                    solut_->setICMoles(i,ICMoles[i]);
                    chemSys_->setDCMoles(waterId,waterMoles);
                    #ifdef DEBUG
                    if (verbose_) {
                        cout << "KineticModel::calculateKineticStep New ICmoles for H is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    }
                    #endif
                }
                if (ICName[i] == "O") {
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Previous IC moles for O is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    #endif
                    ICMoles[i] = waterMoles;
                    solut_->setICMoles(i,ICMoles[i]);
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep New ICmoles for O is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    #endif
               }
            }
            double wmv = chemSys_->getNode()->DC_V0(chemSys_->getDCId("H2O@"),
                                                    chemSys_->getP(),
                                                    chemSys_->getTemperature());
            chemSys_->setGEMPhaseMass(chemSys_->getGEMPhaseId("aq_gen"),waterMass);
            chemSys_->setGEMPhaseVolume(chemSys_->getGEMPhaseId("aq_gen"),wmv/waterMoles);

            if (verbose_) {
                cout << "KineticModel::calculateKineticStep Looping over soluble minerals" << endl;
                cout.flush();
            }
            #ifdef DEBUG
                cout << "KineticModel::calculateKineticStep Looping over soluble minerals" << endl;
                cout << "KineticModel::calculateKineticStep Here is the list of them:" << endl;
                for (int i = 0; i < microPhaseId_.size(); i++) {
                    if (isSoluble(i)) {
                        cout << "KineticModel::calculateKineticStep " << name_[i] << endl;
                    }
                }
                cout.flush();
            #endif
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
                        cout << "KineticModel::calculateKineticStep " << name_[i]
                             << ": Mass dissolved = " << massDissolved << endl;
                        cout.flush();
                        cout << "KineticModel::calculateKineticStep " << name_[i]
                             << ": DC molar mass of "
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
                            cout << "KineticModel::calculateKineticStep " << name_[i]
                                 << ":     Total dissolved "
                                 << ICName[ii] << " = " << ICMoles[ii] << " mol "
                                 << endl;
                        }
                    }
                    if (verbose_) cout.flush();

                } // END OF IF SOLUBLE
            }

            if (verbose_) {
                cout << "KineticModel::calculateKineticStep Looping over thermo minerals" << endl;
                cout.flush();
            }    
            #ifdef DEBUG
                cout << "KineticModel::calculateKineticStep Here is the list of them:" << endl;
                for (int i = 0; i < microPhaseId_.size(); i++) {
                    if (isThermo(i)) {
                        cout << "KineticModel::calculateKineticStep " << name_[i] << endl;
                    }
                }
                cout.flush();
            #endif
            for (int i = 0; i < microPhaseId_.size(); i++) {
                if (isThermo(i)) {
                    for (int ii = 0; ii < ICMoles.size(); ii++) {
                        ICMoles[ii] += ((scaledMass_[i]
                                    / chemSys_->getDCMolarMass(DCId_[i]))
                                    * chemSys_->getDCStoich(DCId_[i],ii));
                        #ifdef DEBUG
                            cout << "KineticModel::calculateKineticStep " << name_[i]
                                 << ":     Total IC " << ICName[ii] << " = "
                                 << ICMoles[ii] << " mol " << endl;
                            cout.flush();
                        #endif
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
                    cout << "KineticModel::calculateKineticStep "
                         << "modifying initial pore solution" << endl;
                    cout.flush();
                }
                #ifdef DEBUG
                    cout << "KineticModel::calculateKineticStep "
                         << "--->Adding " << p->second
                         << " mol/kgw of "
                         << ICName[p->first] << " to initial solution." << endl;
                    cout.flush();
                #endif
                ICMoles[p->first] += (p->second * kgWaterMass);
                p++;
            }

            // @todo BULLARD PLACEHOLDER
            // Still need to implement constant gas phase composition
           
            // @todo BULLARD PLACEHOLDER
            // While we are still using PK model for clinker phases, we
            // must scan for and account for pozzolanic reactive components
            // that will alter the hydration rate of clinker
            // phases, presumably by making a denser hydration shell
            // with slower transport

            for (int i = 0; i < microPhaseId_.size(); i++) {

              microPhaseId = microPhaseId_[i];

              // @todo BULLARD PLACEHOLDER
              // Currently silica fume is the only pozzolanically reactive
              // component
            
              if (chemSys_->getMicroPhaseName(microPhaseId) == "Silica-amorph") {
                  double betarea = k1_[i];
                  double area = betarea * scaledMass_[i];
                  double loi = k2_[i];
                  double sio2 = critDOH_[i];

                  // @note Handle silica fume as a special case here:
                  // Only for silica fume, we let k1 = BET surface area in m2/g,
                  // k2 = LOI in percent by solid mass.
                  // critDOH = silica content in PERCENT BY MASS of silica fume
                  
                  for (int ii = 0; ii < microPhaseId_.size(); ii++) {
                      if (isPK(ii)) {
                          pfk1_ = pfk2_ = pfk3_ = ((betarea/24.0) * (sio2/94.0) * (sio2/94.0));
                      }
                  }
              }
            }

        }   // End of special first-time tasks

        if (hyd_time < leachTime_ && hyd_time < sulfateAttackTime_) { 

          // @todo BULLARD PLACEHOLDER
          // Still need to implement constant gas phase composition
          // Will involve equilibrating gas with aqueous solution
          //
          // First step each iteration is to equilibrate gas phase
          // with the electrolyte, while forbidding anything new
          // from precipitating.
           
          #ifdef DEBUG
             cout << "KineticModel::calculateKineticStep "
                  << "Looping over kinetically controlled phases.  " << endl;
             cout.flush();
          #endif

          vector<double> impurityRelease;
          impurityRelease.clear();
          impurityRelease.resize(chemSys_->getNumMicroImpurities(),0.0);

          // RH factor is the same for all clinker phases
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
          ///    where d is the pore diameter in meters and T is absolute temperature
          
          /// Assume a zero contact angle for now.
          /// @todo revisit the contact angle issue
          
          double critporediam = lattice_->getLargestSaturatedPore(); // in nm
          critporediam *= 1.0e-9;                                    // in m
          double rh = exp(-6.23527e-7/critporediam/T);

          /// Assume a zero contact angle for now.
          /// @todo revisit the contact angle issue
          
          rh = rh > 0.55 ? rh : 0.551;
          rhFactor = pow(((rh - 0.55)/0.45),4.0);
          

          for (int i = 0; i < microPhaseId_.size(); i++) {
            microPhaseId = microPhaseId_[i];
            if (isKinetic(i)) {
                if (initScaledMass_[i] > 0.0) {
                  DOH = (initScaledMass_[i] - scaledMass_[i]) /
                        (initScaledMass_[i]);
                  if (!isPK(i)) {
                      DOH = min(DOH,0.99);  // prevents DOH from prematurely stopping PK calculations
                  }
                } else {
                  throw FloatException("KineticModel","calculateKineticStep",
                             "initScaledMass_ = 0.0");
                }

                cDOH = 1.333 * wcRatio_;
                wcFactor = 1.0;
                if (DOH > cDOH) {
                    wcFactor += ((4.444 * wcRatio_) -
                                 (3.333 * DOH));
                    wcFactor = pow(wcFactor,4.0);
                }

                /*
                wcFactor = 1.0 + (3.333 *
                       (pow(((critDOH_[i] * wcRatio_) - DOH),4.0)));
                */

                arrhenius = exp((activationEnergy_[i]/GASCONSTANT)*((1.0/refT_) - (1.0/T)));

                if (DOH < 1.0 && !doTweak) {
                    
                    // @todo BULLARD PLACEHOLDER
                    // Handle silica fume as a special case here:
                    // Only for silica fume, we let k1 = BET surface area in m2/g,
                    // k2 = LOI in percent by solid mass.
                    // critDOH = silica content in PERCENT BY MASS of silica fume
                    //
                    // This is a TOTAL KLUGE right now from here...

                    if (chemSys_->getMicroPhaseName(microPhaseId) == "Silica-amorph") {

                        // TBD, Dove and Crerar 1990 (Geochim et Cosmochim Acta) 
                        // report k+ (aSiO2)(aH2O)^2 (1 - Q/K) for pure water, where
                        // k+ = 5.6e-9 mol m-2 s-1  in 0.05 m NaCl at 100 C and
                        // an activation enthalpy of about 71 kJ/mol
                        
                        double baserateconst = 5.6e-3 * arrhenius;  // mol m-2 s-1
                                                                    //

                        double ca = chemSys_->getDCConcentration("Ca+2");
                        double kca = 4.0e-7; // mol m-2 s-1 ads. rate const for Ca (guess)
                        double Kca = 10.0; // adsorption equilibrium constant is a guess
                        double na = chemSys_->getDCConcentration("Na+");
                        double kna = 6.35e-7; // mol m-2 s-1 ads. rate const from Dove and Crerar
                        double Kna = 58.3; // adsorption equilibrium constant from Dove and Crerar
                        double k = chemSys_->getDCConcentration("K+");
                        double kk = 5.6e-7; // mol m-2 s-1 ads. rate const from Dove and Crerar
                        double Kk = 46.6; // adsorption equilibrium constant from Dove and Crerar

                        // Langmuir adsorption isotherms assumed to be additive
                       
                        baserateconst += (kca * Kca * ca / (1.0 + (Kca * ca)));
                        baserateconst += (kna * Kna * na / (1.0 + (Kna * na)));
                        baserateconst += (kk * Kk * k / (1.0 + (Kk * k)));
 
                        double ohact = chemSys_->getDCActivity("OH-");
                        double betarea = k1_[i];
                        double area = betarea * scaledMass_[i];
                        double loi = k2_[i];
                        double sio2 = critDOH_[i];
                        arrhenius = exp((activationEnergy_[i]/GASCONSTANT)*((1.0/373.0) - (1.0/T)));
                        double molarmass = chemSys_->getDCMolarMass("Amor-Sl"); // g mol-1
                        vector<int> GEMPhaseIndex;
                        GEMPhaseIndex.clear();
                        GEMPhaseIndex = chemSys_->getMicroPhaseToGEMPhase(microPhaseId);

                        // saturation index , but be sure that there is only one GEM Phase
                       
                        double satindex = solut_->getSI(GEMPhaseIndex[0]);

                        // activity of water
                        double wateractivity = chemSys_->getDCActivity(chemSys_->getDCId("H2O@"));

                        // Make sure sio2 is given in percent by mass of silica fume
                        // This equation basically implements the Dove and Crerar rate equation
                        // for quartz.  Needs to be calibrated for silica fume, but hopefully
                        // the BET area and LOI will help do that.
                       
                        rate = baserateconst * ohact * area * pow(wateractivity,2.0)
                               * (1.0 - (loi/100.0)) * (sio2/100.0) * (1.0 - satindex); 
                                    
                        massDissolved = (rate * timestep) * molarmass;

                        cout << "Silica-amorph: baserateconst = " << baserateconst << " mol/m2/s" << endl;
                        cout << "Silica-amorph: betarea = " << betarea << " m2/g" << endl;
                        cout << "Silica-amorph: scaledMass = " << scaledMass_[i] << " g/(100 g)" << endl;
                        cout << "Silica-amorph: area = " << area << " m2" << endl;
                        cout << "Silica-amorph: loi = " << loi << " %" << endl;
                        cout << "Silica-amorph: sio2 = " << sio2 << " %" << endl;
                        cout << "Silica-amorph: arrhenius = " << arrhenius << endl;
                        cout << "Silica-amorph: satindex = " << satindex << endl;
                        cout << "Silica-amorph: wateractivity = " << wateractivity << endl;
                        cout << "Silica-amorph: rate = " << rate << " mol/s" << endl;
                        cout << "Silica-amorph: massDissolved = " << massDissolved << " g/(100 g)" << endl;
                        cout.flush();

                        scaledMass_[i] -= massDissolved;


                    // @todo BULLARD PLACEHOLDER ... to here
                    
                    } else {
                        cDOH = 1.333 * wcRatio_;
                        wcFactor = 1.0;
                        if (DOH > cDOH) {
                            wcFactor += ((4.444 * wcRatio_) -
                                         (3.333 * DOH));
                            wcFactor = pow(wcFactor,4.0);
                        }

                        // Normal Parrott and Killoh implementation here
                        
                        // @todo BULLARD PLACEHOLDER pfk1_, pfk2_, pfk3_ are
                        // pozzolanic modification factors for the Parrott and Killoh
                        // rate constants k1, k2, and k3.
        
                        if (fabs(n1_[i]) > 0.0) {
                          ngrate = (pfk1_ * k1_[i]/n1_[i]) * (1.0 - DOH)
                                 * pow((-log(1.0 - DOH)),(1.0 - n1_[i]));
                          ngrate *= (blaineFactor_);  // only used for the N+G rate
              
                          if (ngrate < 1.0e-10) ngrate = 1.0e-10;
                        } else {
                          throw FloatException("KineticModel","calculateKineticStep",
                                         "n1_ = 0.0");
                        }
        
                        hsrate = pfk3_ * k3_[i] * pow((1.0 - DOH),n3_[i]);
                        if (hsrate < 1.0e-10) hsrate = 1.0e-10;

                        if (DOH > 0.0) {
                            diffrate = (pfk2_ * k2_[i] * pow((1.0 - DOH),(2.0/3.0))) /
                                       (1.0 - pow((1.0 - DOH),(1.0/3.0)));
                            if (diffrate < 1.0e-10) diffrate = 1.0e-10;
                        } else {
                            diffrate = 1.0e9;
                        }

                        rate = (ngrate < hsrate) ? ngrate : hsrate;
                        if (diffrate < rate) rate = diffrate;
                        rate *= (wcFactor * rhFactor * arrhenius);
                        newDOH = DOH + (rate * timestep);

                        #ifdef DEBUG
                            cout << "KineticModel::calculateKineticStep PK model for " << getName(i)
                                 << ", ngrate = " << ngrate
                                 << ", hsrate = " << hsrate
                                 << ", diffrate = " << diffrate
                                 << ", rhFactor = " << rhFactor
                                 << ", wcFactor = " << wcFactor
                                 << ", timestep " << timestep
                                 << ", oldDOH = " << DOH << ", new DOH = "
                                 << newDOH << endl;
                            cout.flush();
                        #endif

                        /// @note This where we can figure out the volume dissolved
                        /// and link it back to the current volume to see how many
                        /// voxels need to dissolve
                        ///

                        /// @note This all depends on concept of degree of
                        /// hydration as defined by the PK model 

                        /// @todo Make this independent of PK model
                    
                        scaledMass_[i] = initScaledMass_[i] * (1.0 - newDOH);
                        massDissolved = (newDOH - DOH) * initScaledMass_[i];

                        // end it here?
                    }
                       
                    chemSys_->setMicroPhaseMass(microPhaseId,scaledMass_[i]);
                    chemSys_->setMicroPhaseMassDissolved(microPhaseId,massDissolved);

                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Original scaled mass = " << initScaledMass_[i]
                             << " and new scaled mass = "
                             << chemSys_->getMicroPhaseMass(microPhaseId)
                             << " and new volume = "
                             << chemSys_->getMicroPhaseVolume(microPhaseId) << endl;
                        cout.flush();
                    #endif

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

                      /// @todo BULLARD PLACEHOLDER
                      /// Special case for Silica-amorph here to account for SiO2 < 100%
                     
                      double purity = 1.0;
                      if (chemSys_->getMicroPhaseName(microPhaseId) == "Silica-amorph") {
                          purity = critDOH_[i] / 100.0;
                      }

                      ICMoles[ii] += ((purity * massDissolved
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

        #ifdef DEBUG
            if (!doTweak) {
                cout << "KineticModel::calculateKineticStep ICmoles after dissolving:" << endl;
            } else {
                cout << "KineticModel::calculateKineticStep ICmoles after tweaking:" << endl;
            }
            for (int i = 0; i < ICNum; i++) {
              cout << "    " << ICName[i] << ": " << ICMoles[i] << " mol" << endl;
            }
            cout.flush();
        #endif

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
        #ifdef DEBUG
            cout << "KineticModel::calculateKineticStep Aqueous mass = "
                 << chemSys_->getGEMPhaseMass("aq_gen") << endl;
            cout.flush();
        #endif
        for (int i = 0; i < microPhaseId_.size(); i++) {
            microPhaseId = microPhaseId_[i];
            if (isThermo(i)) {
                #ifdef DEBUG
                    cout << "KineticModel::calculateKineticStep Partitioning alkali for phase " << name_[i] << ": ";
                    cout.flush();
                #endif
                phaseList = chemSys_->getMicroPhaseToGEMPhase(microPhaseId);
                // List of all phases making up this microstructure phase
                // Compute the mass of each phase and sum to get phase mass
                phMass = dphMass = 0.0;
                for (int j = 0; j < phaseList.size(); j++) {
                    // Mass will be in g, not kg
                    #ifdef DEBUG
                        cout << endl;
                        cout << "KineticModel::calculateKineticStep phaseName: "
                             << chemSys_->getGEMPhaseName(phaseList[j]) << endl;
                        cout << "KineticModel::calculateKineticStep phaseMass: "
                             << chemSys_->getGEMPhaseMass(phaseList[j]) << endl;
                        cout.flush();
                    #endif
                    dphMass += (double)(chemSys_->getGEMPhaseMass(phaseList[j])
                                      - chemSys_->getPrevGEMPhaseMass(phaseList[j]));
                    phMass += (double)(chemSys_->getGEMPhaseMass(phaseList[j]));
                }
                #ifdef DEBUG
                    cout << "KineticModel::calculateKineticStep Mass is " << phMass << "; mass increment is "
                         << dphMass << endl;
                    cout.flush();
                #endif
                for (int j = 0; j < RdICId_[i].size(); j++) {
                    ICId = RdICId_[i][j];
                    Rd = Rd_[i][j];
                    ICMol = chemSys_->getGEMPhaseStoich("aq_gen",ICId);
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep  Partitioning IC "
                             << chemSys_->getICName(ICId) << " with aqueous moles = ";
                        cout << ICMol << endl;
                        cout << "KineticModel::calculateKineticStep  Partitioning IC " << chemSys_->getICName(ICId);
                        cout << " with actual aqueous moles = ";
                        cout << chemSys_->getGEMPhaseStoich(chemSys_->getGEMPhaseId("aq_gen"),ICId) << endl;
                        cout << "KineticModel::calculateKineticStep moles of 'O' in the solution: ";
                        cout << chemSys_->getGEMPhaseStoich(chemSys_->getGEMPhaseId("aq_gen"),
                                                         chemSys_->getICId("O")) << endl;
                        cout.flush();
                    #endif
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
                            #ifdef DEBUG
                                cout << "KineticModel::calculateKineticStep partitionedMoles for Na: "
                                     << partitionedMoles << endl;
                                cout << "KineticModel::calculateKineticStep actually removed Na: " << partitionedMoles_Na
                                     << " moles" << endl;
                                cout.flush();
                            #endif
                            ICMoles[ICId] = ICMoles[ICId] - partitionedMoles_Na;
                        }
                        if (ICId == chemSys_->getICId("K")) {
                            partitionedMoles_K = partitionedMoles * 
                            chemSys_->getGEMPhaseStoich(chemSys_->getGEMPhaseId("aq_gen"),
                                             chemSys_->getICId("O"));
                            #ifdef DEBUG
                            if (verbose_) {
                                cout << "KineticModel::calculateKineticStep partitionedMoles for K: "
                                     << partitionedMoles << endl;
                                cout << "KineticModel::calculateKineticStep actually removed K: " << partitionedMoles_K
                                     << " moles" << endl;
                            }
                            #endif
                            ICMoles[ICId] = ICMoles[ICId] - partitionedMoles_K;
                        }
                        if (ICMol < 1.0e-16) ICMol = 1.0e-16;
                    }
                }
            }  // END OF IF THERMO
        }  

        if (hyd_time >= leachTime_ && hyd_time < sulfateAttackTime_) {

            ///
            /// This is the block for simulating leaching kinetics
            ///

            if (verbose_) {
                cout << "KineticModel::calculateKineticStep Leaching module at time = "
                     << hyd_time << endl;
                cout.flush();
            }

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
            if (verbose_) {
                cout << "KineticModel::calculateKineticStep Setting target calcium concentration to be: "
                     << targetConc << endl;
            }
    
            if ((chemSys_->getNode())->Get_cDC((long int) CaId) > (targetConc + 0.0001)) {
                Ca2_aq = 1.0e-7 * waterMass * 1.0e-3;
            } else {
                Ca2_aq = targetConc * waterMass * 1.0e-3;
            }
	
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
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << removed_K << "moles of K." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "Na") {
                    ICMoles[ii] = ICMoles[ii] - removed_Na;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << removed_Na << "moles of Na." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "Ca") {
                    ICMoles[ii] = ICMoles[ii] - removed_Ca;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << removed_Ca << "moles of Ca." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "S") {
                    ICMoles[ii] = ICMoles[ii] - removed_S;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << removed_S << "moles of S." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "O") {
                    removed_O = removed_Na + removed_K + 2 * removed_Ca + 4 * removed_S;
                    ICMoles[ii] = ICMoles[ii] - removed_O;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << removed_O << "moles of O." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "H") {
                    removed_H = removed_Na + removed_K + 2 * removed_Ca;
                    ICMoles[ii] = ICMoles[ii] - removed_H;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << removed_H << "moles of H." << endl;
                        cout.flush();
                    #endif
                }
            }

        }   // End of block on leaching kinetics

        for (int ii = 0; ii < ICMoles.size(); ii++) {
            chemSys_->setICMoles(ii,ICMoles[ii]);
            #ifdef DEBUG
                cout << "KineticModel::calculateKineticStep ICmoles of " << ICName[ii]
                     << " is: " << ICMoles[ii] << endl;
                cout.flush();
            #endif
        }

        if (hyd_time >= sulfateAttackTime_ && hyd_time < leachTime_) {
    
            ///
            /// This is the block for simulating sulfate attack kinetics
            ///

            if (verbose_) {
                cout << "KineticModel::calculateKineticStep Sulfate attack module at time = "
                     << hyd_time << endl;
                cout.flush();
            }

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

            double c1; //current concentration
            double c2; //new concentration

            double waterMass = 0.0;
            waterMass = chemSys_->getDCMoles(chemSys_->getDCId("H2O@"))
                        * chemSys_->getDCMolarMass("H2O@");
            
            int KId = chemSys_->getDCId("K+");
            int NaId = chemSys_->getDCId("Na+");
            int CaId = chemSys_->getDCId("Ca+2");
            int CaId1 = chemSys_->getDCId("CaOH+");
            int SId = chemSys_->getDCId("SO4-2");

            /// @todo Better checking for validity of file
            /// and whether it opened
           
            ifstream in("targetConc.dat");
            double concBuff = 0.0;
            in >> concBuff;
            double NaTarget = concBuff;
            in >> concBuff;
            double STarget = concBuff;
            in >> concBuff;
            double CaTarget = concBuff;
            in >> concBuff;
            double Ca1Target = concBuff;
            in >> concBuff;
            double KTarget = concBuff;
            in.close();    
            if (verbose_) {
                cout << "KineticModel::calculateKineticStep NaTarget = " << NaTarget << endl;
                cout << "KineticModel::calculateKineticStep STarget = " << STarget << endl;
                cout << "KineticModel::calculateKineticStep CaTarget = " << CaTarget << endl;
                cout << "KineticModel::calculateKineticStep Ca1Target = " << Ca1Target << endl;
                cout << "KineticModel::calculateKineticStep KTarget = " << KTarget << endl;
                cout.flush();
            }

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
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Added "
                             << add_Na << " moles of Na." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "S") {
                    ICMoles[ii] = ICMoles[ii] + add_S;
                    solutICMoles[ii] = solutICMoles[ii] + add_S;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Added "
                             << add_S << " moles of S." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "K") {
                    ICMoles[ii] = ICMoles[ii] - remove_K;
                    solutICMoles[ii] = solutICMoles[ii] - remove_K;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << remove_K << " moles of K." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "Ca") {
                    ICMoles[ii] = ICMoles[ii] - remove_Ca;
                    solutICMoles[ii] = solutICMoles[ii] - remove_Ca;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Removed "
                             << remove_Ca << " moles of Ca." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "O") {
                    ICMoles[ii] = ICMoles[ii] + add_O - remove_O;
                    solutICMoles[ii] = solutICMoles[ii] + add_O - remove_O;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Added "
                             << (add_O - remove_O) << " moles of O." << endl;
                        cout.flush();
                    #endif
                }
                if (ICName[ii] == "H") {
                    ICMoles[ii] = ICMoles[ii] + add_H - remove_H;
                    solutICMoles[ii] = solutICMoles[ii] + add_H - remove_H;
                    #ifdef DEBUG
                        cout << "KineticModel::calculateKineticStep Added "
                             << (add_H - remove_H) << " moles of H." << endl;
                        cout.flush();
                    #endif
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

    #ifdef DEBUG
        cout << "KineticModel::calculatePhaseChange Microphase "
             << chemSys_->getMicroPhaseName(microPhaseId) << " contains phases: "
             << endl;
        for (int i = 0; i < GEMPhaseIndex.size(); i++) {
            cout << "KineticModel::calculatePhaseChange            "
                 << chemSys_->getGEMPhaseName(i) << endl;
        }    
        cout.flush();
    #endif

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
      #ifdef DEBUG
          cout << "KineticModel::calculatePhaseChange Phase "
               << chemSys_->getGEMPhaseName(GEMPhaseIndex[i]) 
               << " contains DC members: " << endl;
      #endif
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

    #ifdef DEBUG
        cout << "KineticModel::setKineMicDCmoles" << endl;
        cout.flush();
    #endif

    try {
        int waterId = chemSys_->getDCId("H2O@");
        double waterMoles = chemSys_->getDCMoles(waterId);
        double waterMolarMass = chemSys_->getDCMolarMass(waterId);
        double waterMass = waterMoles * waterMolarMass;
        #ifdef DEBUG
            cout << "KineticModel::setKineticDCmoles        "
                 << chemSys_->getDCName(waterId) << ": Mass = " << waterMass
                 << ", Molar mass = " << waterMolarMass << endl;
            cout.flush();
        #endif
        for (int i = 0; i < microPhaseId_.size(); i++) {
            if (isKinetic(i)) {
                if (chemSys_->getDCMolarMass(DCId_[i]) <= 0.0) {
                    throw FloatException("KineticModel","setKineticDCmoles",
                                         "Divide by zero error");
                }
                #ifdef DEBUG
                    cout << "KineticModel::setKineticDCmoles        Clinker phase "
                         << name_[i] << ": Mass = " << scaledMass_[i]
                         << ", Molar mass = " << chemSys_->getDCMolarMass(DCId_[i])
                         << endl;
                #endif
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
