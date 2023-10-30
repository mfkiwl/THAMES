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
                
                parsePhase(doc, cur, numEntry);
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

            makeModel(doc, cur, kineticData);

            /// STOPPED HERE
        }
        cur = cur->next;
    }

    if (kineticfound) {

        if (kineticData.type == "parrotkilloh") {
            
            #ifdef DEBUG
                cout << "KineticControllerl::parsePhase PK " << kineticData.name
                     << ", id = "
                     << kineticData.microPhaseId << endl;
                cout << "KineticController::parsePhase   (k1,k2,k3) =  " << kineticData.k1
                     << "," << kineticData.k2
                     << "," << kineticData.k3 << endl;
                cout << "KineticController::parsePhase   (n1,n3) =  " << kineticData.n1 << ","
                     << kineticData.n3 << endl;
                cout << "KineticModel::parsePhase   Ea = " << kineticData.Ea << endl;
                cout.flush();
            #endif
            isParrotKilloh = true;
            ispozz = false;
            istherm = false;
            issol = false;

        } else if (kineticData.type == "pozzolanic") {
            #ifdef DEBUG
                cout << "KineticModel::parsePhase pozzolanic " << kineticData.name
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
            isParrotKilloh = false;
            ispozz = true;
            istherm = false;
            issol = false;


        } else if (kineticData.type == "soluble") {

            iskin = false;
            ispozz = false;
            istherm = true;
            issol = true;

        } else {

            iskin = false;
            ispozz = false;
            istherm = true;
            issol = false;
        }
    } else {
        isParrotKilloh = ispozz = istherm = issol = false;
    }
    
    isParrotKilloh_.push_back(isParrotKilloh);
    isPozzolanic_.push_back(ispozz);
    isThermo_.push_back(istherm);
    isSoluble_.push_back(issol);

    // @todo BULLARD PLACEHOLDER
    // Special kluges here for silica fume, which we have to call Silica-amorph
    // for now to make it compatible with the user interface as of 2022 Dec 29
    //
    // @todo Find a way to make this general to all pozzolans
   
    string sfume("SFUME");
    string siamorph("Silica-amorph");
    if (kineticData.name == sfume || kineticData.name == siamorph) {
        kineticData.name = siamorph;
        isParrotKilloh_.push_back(false);
    } else {
        isParrotKilloh_.push_back(iskin);
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

void KineticController::makeModel (xmlDocPtr doc,
                                   xmlNodePtr cur,
                                   int &numEntry,
                                   KineticData &kineticData)
{
    KineticModel *km = NULL;

    if (kineticData.type == "parrotkilloh") {
        km = new ParrotKillohModel(kineticData);
    } else if (kineticData.type == "pozzolanic") {
        km = new PozzolanicModel(kineticData);
    } else if (kineticData.type == "ordinary") {
        km = new OrdinaryModel(kineticData);
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

