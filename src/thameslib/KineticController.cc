/**
@file  KineticController.cc
@brief Method definitions for the KineticController class.

*/
#include "KineticController.h"

KineticController::KineticController ()
{
    temperature_ = 293.15;  // default temperature (K)
    refT_ = 293.15;         // default temperature (K)

    ///
    /// Clear out the vectors so they can be populated with values from the
    ///

    numPhases_ = 0;
    chemSys_ = NULL;
    solut_ = NULL;
    lattice_ = NULL;
    phaseKineticModel_.clear();
    name_.clear();
    isKinetic_.clear();
    
    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10;

    verbose_ = warning_ = false;

    return;
}

KineticController::KineticController (ChemicalSystem *cs,
                                      Solution *solut,
                                      Lattice *lattice,
                                      const string &fileName,
                                      const bool verbose,
                                      const bool warning)
:chemSys_(cs),solut_(solut),lattice_(lattice)
{
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
        cout << "KineticController::KineticController Finished reading chemistry.xml file" << endl;
        for (int i = 0; i < microPhaseId_.size(); ++i) {
            microPhaseId = microPhaseId_[i];
            if (isKinetic(i)) {
                cout << "KineticController::KineticController kinetic phase " << microPhaseId << endl;
                cout << "KineticController::KineticController     name = " << chemSys_->getMicroPhaseName(microPhaseId)
                     << endl;
            }
        }
        cout.flush();
    #endif

    // The ChemicalSystem and Lattice objects were constructed before we get
    // here, so now we know nearly everything we need to know about the
    // microstructure.  The only thing left to do is to populate the scaled
    // phase masses and read the w/s ratio, which the kinetic model needs
    
    getPhaseMasses();

    return;
}

void KineticController::parseDoc (const string &docName)
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
            if ((!xmlStrcmp(cur->name, (const xmlChar *)"temperature"))) {
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

void KineticController::parsePhase (xmlDocPtr doc,
                                    xmlNodePtr cur,
                                    int &numEntry
                                    KineticData &kineticData)
{
    xmlChar *key;
    int proposedgemphaseid,proposedDCid;
    int testgemid,testdcid;
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
            isKinetic_[isKinetic_.size()-1] = true;
            kineticData.GEMPhaseId =
                chemSys_->getMicroPhaseToGEMPhase(kineticData.microPhaseId,0);
            kineticData.DCId =
                chemSys_->getMicroPhaseDCMembers(kineticData.microPhaseId,0);

            ///
            /// Kinetic data are grouped together,
            /// so there is a method written just for parsing that grouping
            ///

            parseKineticData(doc, cur, kineticData);
            kineticData.temperature = temperature_;
            kineticData.reftemperature = refT_;

            makeModel(doc, cur, kineticData);

        }
        cur = cur->next;
    }


    return;
}

void KineticModel::parseKineticData (xmlDocPtr doc,
                                     xmlNodePtr cur,
                                     KineticData &kineticData)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"type"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            kineticData.type = st;
            xmlFree(key);
        }
        // Parrot-Killoh k1 parameter
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k1"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.k1,st);
            xmlFree(key);
        }
        // Parrot-Killoh k2 parameter
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k2"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.k2,st);
            xmlFree(key);
        }
        // Parrot-Killoh k3 parameter
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k3"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.k3,st);
            xmlFree(key);
        }
        // Parrot-Killoh n1 parameter
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"n1"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.n1,st);
            xmlFree(key);
        }
        // Parrot-Killoh n3 parameter
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"n3"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.n3,st);
            xmlFree(key);
        }
        // Parrot-Killoh critical DOH parameter
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"critdoh"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.critDOH,st);
            xmlFree(key);
        }
        // Generic flux-like rate constant for dissolution or precipitation
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"rateconst"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.rateconst,st);
            xmlFree(key);
        }
        // Exponent on  the saturation index in the rate equation
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"siexp"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.rateconst,st);
            xmlFree(key);
        }
        // Exponent on  the driving force term in the rate equation
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"dfexp"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.rateconst,st);
            xmlFree(key);
        }
        // Exponent on  the hydroxy ion activity in the rate equation
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"ohexp"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.rateconst,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Ea"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kineticData.Ea,st);
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

void KineticController::makeModel (xmlDocPtr doc,
                                   xmlNodePtr cur,
                                   int &numEntry,
                                   KineticData &kineticData)
{
    KineticModel *km = NULL;

    if (kineticData.type == "parrotkilloh") {
        // Read remaining Parrot and Killoh model parameters
        km = new ParrotKillohModel(chemSys_,solut_,lattice_,kineticData,verbose_,warning_);
    } else if (kineticData.type == "pozzolanic") {
        // Read remaining pozzolanic model parameters
        km = new PozzolanicModel(chemSys_,solut_,lattice_,kineticData,verbose_,warning_);
    }

    phaseKineticModel_.push_back(km);

    return;
}

void KineticController::calculateKineticStep (const double timestep,
                                              const double temperature,
                                              bool isFirst)
{
    cout << "Whoop!!" << endl;
}

