/**
@file
@brief Method definitions for the KineticModel class.

*/
#include "KineticModel.h"

KineticModel::KineticModel ()
{
    ///
    /// Default value for w/c ratio in PK model is 0.45
    ///

    wcratio_ = 0.45;

    ///
    /// Default value for Blaine fineness in PK model is 385 m<sup>2</sup>/kg
    ///

    blaine_ = 385.0;
    refblaine_ = 385.0;     // reference Blaine fineness (m2/kg)

    ///
    /// Default temperaturein the PK model is 20 C (or 293 K)
    ///

    temperature_ = 293.15;  // default temperature (K)
    refT_ = 293.15;         // default temperature (K)

    ///
    /// Clinker phases have special status in the PK kinetic model, and occupy
    /// phase id numbers 0, 1, 2, and 3 for alite, belite, aluminate, and ferrite,
    /// respectively
    ///

    modelc3sid_ = 0;
    modelc2sid_ = 1;
    modelc3aid_ = 2;
    modelc4afid_ = 3;

    ///
    /// Clear out the vectors so they can be populated with values from the
    /// XML input file
    ///

    name_.clear();
    kineticphase_.clear();
    thermophase_.clear();
    solublephase_.clear();
    chemsysDCid_.clear();
    chemsysphaseid_.clear();
    RdICid_.clear();
    k1_.clear();
    k2_.clear();
    k3_.clear();
    n1_.clear();
    n3_.clear();
    Rd_.clear();
    Ea_.clear();
    scaledmass_.clear();
    initscaledmass_.clear();
    critDOH_.clear();
    doh_.clear();
    Na_target_.clear();
    K_target_.clear();
    Mg_target_.clear();
    SO4_target_.clear();
    
    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sattack_time_ = 1.0e10;
    leach_time_ = 1.0e10;

    return;
}

KineticModel::KineticModel (ChemicalSystem *cs,
                            Solution *solut,
                            Lattice *lattice,
                            const string &fname,
                            const bool verbose,
                            const bool warning)
:chemsys_(cs),solut_(solut),lattice_(lattice)
{
    const string PHASENUM = "phasenum";
    const string BEGINPHASE = "<phase>";
    const string ENDPHASE = "</phase>";
    const string PHASENAME = "name";
    const string TYPE = "type";
    const string GEMNAME = "gemname";
    const string DCNAME = "dcname";
    const string C3S = "C3S";
    const string C2S = "C2S";
    const string C3A = "C3A";
    const string C4AF = "C4AF";
    const string C4AFGEMNAME = "Ferrite";
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

    ///
    /// Default value for Blaine fineness in PK model is 385 m<sup>2</sup>/kg
    ///

    blaine_ = 385.0;
    refblaine_ = 385.0;
    blainefactor_ = 1.0;    // ratio of blaine_ to refblaine_

    ///
    /// Default value for w/c ratio in PK model is 0.45
    ///

    wcratio_ = 0.45;

    ///
    /// Default temperaturein the PK model is 20 C (or 293 K)
    ///

    temperature_ = 293.15;
    refT_ = 293.15;

    ///
    /// Clinker phases have special status in the PK kinetic model, and occupy
    /// phase id numbers 0, 1, 2, and 3 for alite, belite, aluminate, and ferrite,
    /// respectively
    ///

    modelc3sid_ = 0;
    modelc2sid_ = 1;
    modelc3aid_ = 2;
    modelc4afid_ = 3;

    ///
    /// Clear out the vectors so they can be populated with values from the
    /// XML input file
    ///

    name_.clear();
    kineticphase_.clear();
    thermophase_.clear();
    solublephase_.clear();
    chemsysDCid_.clear();
    chemsysphaseid_.clear();
    RdICid_.clear();
    k1_.clear();
    k2_.clear();
    k3_.clear();
    n1_.clear();
    n3_.clear();
    Rd_.clear();
    Ea_.clear();
    scaledmass_.clear();
    initscaledmass_.clear();
    critDOH_.clear();
    doh_.clear();
    Na_target_.clear();
    K_target_.clear();
    Mg_target_.clear();
    SO4_target_.clear();
    
    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sattack_time_ = 1.0e10;
    leach_time_ = 1.0e10;
    
    ///
    /// Open the input XML file for kinetic data and parse it
    ///

    string xmlext = ".xml";
    size_t foundxml;
    foundxml = fname.find(xmlext);
    try {
      if (foundxml != string::npos) {
          if (verbose_) cout << "KineticModel data file is an XML file" <<endl;
          parseDoc(fname);
      } else {
          throw FileException("KineticModel","KineticModel",fname,
                            "NOT in XML format");
      }
    }
    catch (FileException fex) {
      fex.printException();
      exit(1);
    }
    return;
}

void KineticModel::parseDoc (const string &docname)
{
    int numentry = -1;
    int testgemid;

    ///
    /// The kineticdata structure is used to temporarily hold parsed data
    /// for a given phase before the data are loaded permanently into class members.
    ///

    KineticData kineticdata;

    ///
    /// This method uses the libxml library, so it needs to be added and linked
    /// at compile time.
    ///

    xmlDocPtr doc;
    xmlChar *key;
    xmlNodePtr cur;

    cout.flush();
    doc = xmlParseFile(docname.c_str());

    ///
    /// Check if the xml file is valid and parse it if so.
    ///

    try {
        string rxcsd = xsd_files_path;
        rxcsd+="/chemistry.xsd";
        if (verbose_) cout << "Chemistry xsd file is at " << rxcsd << endl;
        if(!is_xml_valid(doc,rxcsd.c_str())) {
            throw FileException("KineticModel","KineticModel",docname,
                                "xml NOT VALID");
        }

        if (doc == NULL ) {
            throw FileException("KineticModel","KineticModel",docname,
                            "xml NOT parsed successfully");
        }

        cur = xmlDocGetRootElement(doc);

        if (cur == NULL) {
            xmlFreeDoc(doc);
            throw FileException("KineticModel","KineticModel",docname,
                            "xml document is empty");
        }

        cur = cur->xmlChildrenNode;
        int testnumentries;
        while (cur != NULL) {
            if (verbose_) cout << "Key name = " << cur->name << endl;
            if ((!xmlStrcmp(cur->name, (const xmlChar *)"blaine"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(blaine_,st);
                if (verbose_) cout << "Blaine value in interface xml file is "
                                   << blaine_ << endl;
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"refblaine"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(refblaine_,st);
                if (verbose_) cout << "Reference Blaine value in interface xml file is "
                                   << refblaine_ << endl;
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"wcRatio"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(wcratio_,st);
                if (verbose_) cout << "wcRatio value in interface xml file is "
                                   << wcratio_ << endl;
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"temperature"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(temperature_,st);
                if (verbose_) cout << "Temperature in interface xml file is "
                                   << temperature_ << endl;
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"reftemperature"))) {
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string st((char *)key);
                from_string(refT_,st);
                if (verbose_) cout << "Reference temperature value in interface xml file is "
                                   << refT_ << endl;
                xmlFree(key);
            } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"phase"))) {
                if (verbose_) cout << "Preparing to parse a phase from interface xml file..."
                                   << endl;

                /// Each phase is a more complicated grouping of data that
                /// has a separate method for parsing.
                ///

                parsePhase(doc, cur, numentry, kineticdata);
            }
            cur = cur->next;
        }

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
                               int &numentry,
                               KineticData &kdata)
{
    xmlChar *key;
    int proposedgemphaseid,proposedgemdcid;
    int testmicid,testgemid,testdcid;
    string testname;
    bool kineticfound = false;


    if (verbose_) cout << "In function KineticModel::parsePhase" << endl;
    kdata.rdid.clear();
    kdata.rdval.clear();

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if (verbose_) cout << "    Key name = " << cur->name << endl;
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"thamesname"))) {
            if (verbose_) cout << "    Trying to get phase name" << endl;
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string testname((char *)key);
            kdata.name = testname;
            kdata.micid = chemsys_->getMicid(testname);
            if (verbose_) cout << "    Phase name = " << kdata.name << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"kinetic_data"))) {
            if (verbose_) cout << "    Parsing kinetic data..." << endl;
            numentry += 1;
            kineticfound = true;
            kdata.gemphaseid = chemsys_->getMic2phase(kdata.micid,0);
            kdata.gemdcid = chemsys_->getMic2DC(kdata.micid,0);

            ///
            /// The kinetic data is another complicated grouping of data,
            /// so there is a method written just for parsing that grouping
            ///

            parseKineticData(doc, cur, kdata);
            if (verbose_) cout << "    Done parsing kinetic data." << endl;
        }
        cur = cur->next;
    }

    if (kineticfound) {

        if (kdata.type == "kinetic") {
            kineticphase_.push_back(numentry);
            chemsys_->setKineticphase(kdata.micid);
            chemsys_->setMic2kinetic(kdata.micid,numentry);
            if (verbose_) {
                cout << "    Setting Mic2kinetic(" << kdata.micid << ","
                     << numentry << ")... ";
                cout.flush();
                cout << "Done!" << endl;
                cout.flush();
                cout << "    Setting k1_[" << numentry << "] = "
                     << kdata.k1 << endl;
                cout.flush();
            }
            k1_.push_back(kdata.k1);
            if (verbose_) {
                cout << "    Okay, k1_[" << numentry << "] = "
                     << k1_[numentry] << endl;
                cout << "    Setting k2_[" << numentry << "] = "
                     << kdata.k2 << endl;
                cout.flush();
            }
            k2_.push_back(kdata.k2);
            if (verbose_) {
                cout << "    Okay, k2_[" << numentry << "] = "
                     << k2_[numentry] << endl;
                cout << "    Setting k3_[" << numentry << "] = "
                     << kdata.k3 << endl;
                cout.flush();
            }
            k3_.push_back(kdata.k3);
            if (verbose_) {
                cout << "    Okay, k3_[" << numentry << "] = "
                     << k3_[numentry] << endl;
                cout << "    Setting n1_[" << numentry << "] = "
                     << kdata.n1 << endl;
                cout.flush();
            }
            n1_.push_back(kdata.n1);
            if (verbose_) {
                cout << "    Okay, n1_[" << numentry << "] = "
                     << n1_[numentry] << endl;
                cout << "    Setting n3_[" << numentry << "] = "
                     << kdata.n3 << endl;
                cout.flush();
            }
            n3_.push_back(kdata.n3);
            if (verbose_) {
                cout << "    Okay, n3_[" << numentry << "] = "
                     << n3_[numentry] << endl;
                cout << "    Setting Ea_[" << numentry << "] = "
                     << kdata.Ea << endl;
                cout.flush();
            }
            Ea_.push_back(kdata.Ea);

        } else if (kdata.type == "soluble") {

            if (verbose_) {
                cout << "This is a soluble phase, adding "
                     << numentry << " to the list ";
                cout << "of soluble phases" << endl;
            }
            solublephase_.push_back(numentry);
            if (verbose_) {
                cout << "    Setting Mic2kinetic(" << kdata.micid << ","
                     << numentry << ")... ";
                cout.flush();
            }

            ///
            /// There is some discrepancy about whether an instantly
            /// dissolving phase should be considered as kinetically
            /// controlled or thermodynamically controlled

           // chemsys_->setKineticphase(kdata.micid);
           // chemsys_->setMic2kinetic(kdata.micid,numentry);
            chemsys_->setThermophase(kdata.micid);
            chemsys_->setMic2thermo(kdata.micid,numentry);

        } else {

            if (verbose_) {
                cout << "This is a thermo phase, adding "
                     << numentry << " to the list ";
                cout << "of thermo phases" << endl;
            }
            thermophase_.push_back(numentry);
            chemsys_->setThermophase(kdata.micid);
            chemsys_->setMic2thermo(kdata.micid,numentry);
        }
    
        if (verbose_) {
            cout << "    Setting initscaledmass_[" << numentry << "] = "
                 << kdata.scaledmass << endl;
            cout.flush();
        }
        initscaledmass_.push_back(kdata.scaledmass);
        if (verbose_) {
            cout << "    Setting scaledmass_[" << numentry << "] = "
                 << kdata.scaledmass << endl;
            cout.flush();
        }
        scaledmass_.push_back(kdata.scaledmass);
        if (verbose_) {
            cout << "    Setting name_[" << numentry << "] = "
                 << kdata.name << endl;
            cout.flush();
        }
        name_.push_back(kdata.name);
        if (verbose_) {
            cout << "    Setting Micphasemass, micid = "
                 << kdata.micid << ", mass = " << kdata.scaledmass << endl;
        }
        chemsys_->setMicphasemass(kdata.micid,kdata.scaledmass);
        if (verbose_) {
            cout << "    Setting Micphasemassdissolved, micid = "
                 << kdata.micid << endl;
        }
        chemsys_->setMicphasemassdissolved(kdata.micid,0.0);
        if (verbose_) {
            cout << "    Setting chemsysDCid_[" << numentry << "] = "
                 << kdata.gemdcid << endl;
            cout.flush();
        }
        chemsysDCid_.push_back(kdata.gemdcid);
        if (verbose_) {
            cout << "    Setting chemsysphaseid_[" << numentry << "] = "
                 << kdata.gemphaseid << endl;
            cout.flush();
        }
        chemsysphaseid_.push_back(kdata.gemphaseid);
        if (verbose_) {
            cout << "    Setting critDOH_[" << numentry << "] = "
                 << kdata.critdoh << endl;
            cout.flush();
        }

        critDOH_.push_back(kdata.critdoh);
        RdICid_.push_back(kdata.rdid);
        Rd_.push_back(kdata.rdval);
    }

    return;
}

void KineticModel::parseKineticData (xmlDocPtr doc,
                                     xmlNodePtr cur,
                                     KineticData &kdata)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"type"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            kdata.type = st;
            if (verbose_) cout << "        type = " << kdata.type << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"scaledmass"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.scaledmass,st);
            if (verbose_) cout << "        scaledmass = " << kdata.scaledmass << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k1"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.k1,st);
            if (verbose_) cout << "        k1 = " << kdata.k1 << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k2"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.k2,st);
            if (verbose_) cout << "        k2 = " << kdata.k2 << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k3"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.k3,st);
            if (verbose_) cout << "        k3 = " << kdata.k3 << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"n1"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.n1,st);
            if (verbose_) cout << "        n1 = " << kdata.n1 << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"n3"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.n3,st);
            if (verbose_) cout << "        n3 = " << kdata.n3 << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Ea"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.Ea,st);
            if (verbose_) cout << "        Ea = " << kdata.Ea << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"critdoh"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(kdata.critdoh,st);
            if (verbose_) cout << "        critdoh = " << kdata.critdoh << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Rd"))) {
            if (verbose_) cout << "    Parsing Rd data..." << endl;

            ///
            /// The data about partitioning of impurities among the clinker
            /// phases are grouped within a complex field in the input XML
            /// file, so we have a special method to parse it.
            ///

            parseRdData(doc, cur, kdata);
            if (verbose_) cout << "    Done parsing Rd data." << endl;
        }
        cur = cur->next;
    }

    return;
}

void KineticModel::parseRdData(xmlDocPtr doc,
                               xmlNodePtr cur,
                               KineticData &kdata) 
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;
    int rdid;
    double rdval;

    while (cur != NULL) {
        if (verbose_) cout << "    Key name = " << cur->name << endl;
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Rdelement"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            rdid = chemsys_->getICid(st);
            kdata.rdid.push_back(rdid);
            if (verbose_) cout << "            rdid = " << rdid << endl;
            xmlFree(key);
        }

        if ((!xmlStrcmp(cur->name, (const xmlChar *)"Rdvalue"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(rdval,st);
            kdata.rdval.push_back(rdval);
            if (verbose_) cout << "            rdval = " << rdval << endl;
            xmlFree(key);
        }
        cur = cur->next;
    }
}

void KineticModel::calculateKineticStep (const double timestep,
                                         const double temperature,
                                         bool isfirst)
{
    ///
    /// Initialize local variables
    ///

    double T = temperature;
    double wcfactor = 1.0;
    double DOH = 0.0;
    double arrhenius = 1.0;


    double ngrate = 1.0e-10;          // Nucleation and growth rate
    double hsrate = 1.0e-10;          // Hydration shell rate
    double diffrate = 1.0e-10;        // Diffusion rate
    double rate = 1.0e-10;            // Selected rate

    double newDOH = 0.0;              // Updated value of doh
    double massdissolved = 0.0;

    static double hyd_time = 0.0;
    hyd_time = hyd_time + timestep;
    if (verbose_) cout << "hyd_time = " << hyd_time << endl;

    try {
        static int conc_index = 0;     
        int icnum = chemsys_->getICnum();
        int dcnum = chemsys_->getDCnum();
        int phasenum = chemsys_->getPhasenum();
        int kpid,spid,tpid,dcid,icid;
        double molarmass,Rd;
        vector<double> icmoles,solut_icmoles,dcmoles,phasemoles;
        icmoles.clear();
        icmoles.resize(icnum,0.0);
        solut_icmoles.clear();
        solut_icmoles.resize(icnum,0.0);
        dcmoles.clear();
        dcmoles.resize(dcnum,0.0);
        phasemoles.clear();
        phasemoles.resize(phasenum,0.0);
        string icn;

        vector<string> icname;
        icname.clear();
        icname.resize(icnum," ");
        for (int i = 0; i < icnum; i++) {
          icmoles[i] = chemsys_->getICmoles(i);
          icname[i] = chemsys_->getICname(i);
        }
        for (int i = 0; i < dcnum; i++) {
          dcmoles[i] = chemsys_->getDCmoles(i);
        }
        for (int i = 0; i < phasenum; i++) {
          phasemoles[i] = chemsys_->getPhasemoles(i);
        }

        solut_icmoles = chemsys_->getSolution();

        if (hyd_time < leach_time_ && hyd_time < sattack_time_) { 
          if (verbose_) {
             cout << "Looping over clinker minerals.  " << endl;
             cout.flush();
          }

          vector<double> impurityrelease;
          impurityrelease.clear();
          impurityrelease.resize(chemsys_->getMicimpuritynum(),0.0);

          for (int i = 0; i < kineticphase_.size(); i++) {
            kpid = kineticphase_[i];
            wcfactor = 1.0 + (3.333 *
                   (pow(((critDOH_[kpid] * wcratio_) - DOH),4.0)));

            arrhenius = exp((Ea_[kpid]/GASCONSTANT)*((1.0/refT_) - (1.0/T)));
            if (verbose_) {
                cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ" << endl;
                cout << "Calculating Arrhenius effect:" << endl;
                cout << "    GASCONSTANT = " << GASCONSTANT << endl;
                cout << "    Ea_[" << kpid << "] = " << Ea_[kpid] << endl;
                cout << "    refT_ = " << refT_ << endl;
                cout << "    T = " << T << endl;
                cout << "    arrhenius factor = " << arrhenius << endl;
                cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ" << endl;
            }

            if (initscaledmass_[kpid] > 0.0) {
              DOH = (initscaledmass_[kpid] - scaledmass_[kpid]) /
                    (initscaledmass_[kpid]);
            } else {
              throw FloatException("KineticModel","calculateKineticStep",
                             "initscaledmass_ = 0.0");
            }

            if (DOH < 1.0) {
        
              if (DOH < 1.0) {
                if (fabs(n1_[kpid]) > 0.0) {
                  ngrate = (k1_[kpid]/n1_[kpid]) * (1.0 - DOH)
                         * pow((-log(1.0 - DOH)),(1.0 - n1_[kpid]));
                  ngrate *= (blainefactor_);  // only used for the N+G rate
              
                  if (ngrate < 1.0e-10) ngrate = 1.0e-10;
                } else {
                  throw FloatException("KineticModel","calculateKineticStep",
                                 "n1_ = 0.0");
                }
              } else {
                throw DataException("KineticModel","calculateKineticStep",
                            "DOH >= 1.0");
              }
        
            if (DOH < 1.0) {
              hsrate = k3_[kpid] * pow((1.0 - DOH),n3_[kpid]);
              if (hsrate < 1.0e-10) hsrate = 1.0e-10;
            } else {
               throw DataException("KinetiModel","calculateKineticStep",
                            "DOH >= 1.0");
            }

            if (DOH < 1.0) {
              diffrate = (k2_[kpid] * pow((1.0 - DOH),(2.0/3.0))) /
                         (1.0 - pow((1.0 - DOH),(1.0/3.0)));
              if (diffrate < 1.0e-10) diffrate = 1.0e-10;
            } else {
              throw DataException("KinetiModel","calculateKineticStep",
                            "DOH >= 1.0");
            }


            rate = (ngrate < hsrate) ? ngrate : hsrate;
            if (diffrate < rate) rate = diffrate;
            rate *= (wcfactor * arrhenius);
            newDOH = DOH + (rate * timestep);
            scaledmass_[kpid] = initscaledmass_[kpid] * (1.0 - newDOH);
            massdissolved = (newDOH - DOH) * initscaledmass_[kpid];
            chemsys_->setMicphasemass(chemsys_->getKinetic2mic(kpid),
                              scaledmass_[kpid]);
            chemsys_->setMicphasemassdissolved(chemsys_->getKinetic2mic(kpid),
                                        massdissolved);

            impurityrelease[0] = (massdissolved * chemsys_->getK2o(chemsys_->getKinetic2mic(kpid)));
            impurityrelease[1] = (massdissolved * chemsys_->getNa2o(chemsys_->getKinetic2mic(kpid)));
            impurityrelease[2] = (massdissolved * chemsys_->getMgo(chemsys_->getKinetic2mic(kpid)));
            impurityrelease[3] = (massdissolved * chemsys_->getSo3(chemsys_->getKinetic2mic(kpid)));

            for (int ii = 0; ii < icmoles.size(); ii++) {
              icmoles[ii] += ((massdissolved
                          / chemsys_->getDCmolarmass(chemsysDCid_[kpid]))
                          * chemsys_->getDCstoich(chemsysDCid_[kpid],ii));
              if (icname[ii] == "O") {
                // Dissolved K2O in this clinker phase
                if (chemsys_->isIC("K")) {
                    icn = "K";
                    molarmass = 2.0 * chemsys_->getICmolarmass(icn);
                    icn = "O";
                    molarmass += chemsys_->getICmolarmass(icn);
                    icmoles[ii] += (impurityrelease[0]/molarmass);
                }
                // Dissolved Na2O in this clinker phase
                if (chemsys_->isIC("Na")) {
                    icn = "Na";
                    molarmass = 2.0 * chemsys_->getICmolarmass(icn);
                    icn = "O";
                    molarmass += chemsys_->getICmolarmass(icn);
                    icmoles[ii] += (impurityrelease[1]/molarmass);
                }
                // Dissolved MgO in this clinker phase
                if (chemsys_->isIC("Mg")) {
                    icn = "Mg";
                    molarmass = chemsys_->getICmolarmass(icn);
                    icn = "O";
                    molarmass += chemsys_->getICmolarmass(icn);
                    icmoles[ii] += (impurityrelease[2]/molarmass);
                }
                // Dissolved SO3  in this clinker phase
                if (chemsys_->isIC("S")) {
                    icn = "S";
                    molarmass = chemsys_->getICmolarmass(icn);
                    icn = "O";
                    molarmass += (3.0 * chemsys_->getICmolarmass(icn));
                    icmoles[ii] += (3.0 * (impurityrelease[3]/molarmass));
                }
              } else if (icname[ii] == "S") {
                // Dissolved SO3  in this clinker phase
                icn = "S";
                molarmass = chemsys_->getICmolarmass(icn);
                icn = "O";
                molarmass += (3.0 * chemsys_->getICmolarmass(icn));
                icmoles[ii] += (impurityrelease[3]/molarmass);
              } else if (icname[ii] == "K") {
                // Dissolved K2O in this clinker phase
                icn = "K";
                molarmass = 2.0 * chemsys_->getICmolarmass(icn);
                icn = "O";
                molarmass += chemsys_->getICmolarmass(icn);
                icmoles[ii] += (2.0 * (impurityrelease[0]/molarmass));
              } else if (icname[ii] == "Na") {
                // Dissolved Na2O in this clinker phase
                icn = "Na";
                molarmass = 2.0 * chemsys_->getICmolarmass(icn);
                icn = "O";
                molarmass += chemsys_->getICmolarmass(icn);
                icmoles[ii] += (2.0 * (impurityrelease[1]/molarmass));
              } else if (icname[ii] == "Mg") {
                // Dissolved MgO in this clinker phase
                icn = "Mg";
                molarmass = chemsys_->getICmolarmass(icn);
                icn = "O";
                molarmass += chemsys_->getICmolarmass(icn);
                icmoles[ii] += (impurityrelease[2]/molarmass);
              }
    
            }

            } else {
              throw DataException("KineticModel","calculateKineticStep",
                                "DOH >= 1.0");
            }   
          }
        }	

        if (isfirst) {  // Beginning of special first-time setup tasks
            
            // Set the proper amount of water for the total solid mass
            // and the water-solid ratio
            
            // Determine total intial solid mass
            double totsolidmass = 0.0;
            vector<double> initsolidmass = getInitscaledmass();
            for (int i = 0; i < initsolidmass.size(); i++) {
                totsolidmass += initsolidmass[i];
            }
            double watermass = totsolidmass * getWcratio();
            int waterid = chemsys_->getDCid("H2O@");
            double watermolarmass = chemsys_->getDCmolarmass(waterid);
            double watermoles = watermass / watermolarmass;
            if (verbose_) {
                cout << "*** Initial solid mass = " << totsolidmass << endl;
                cout << "*** w/s ratio = " << getWcratio() << endl;
                cout << "*** Initial water mass = " << watermass << endl;
                cout << "*** Initial water moles = " << watermoles << endl;
                cout << "***" << endl;
            }

            if (watermass <= 0.0) {
                throw FloatException("Controller","calculateState","Divide by zero error");
            }
            for (int i = 0; i < icmoles.size(); i++) {
                if (icname[i] == "H") {
                    if (verbose_) cout << "previous IC moles for H is: " << icmoles[i] << endl;
                    icmoles[i] = (2.0 * watermoles);
                    solut_->setICmoles(i,(2.0 * watermoles));
                    if (verbose_) cout << "new ICmoles for H is: " << icmoles[i] << endl;
                }
                if (icname[i] == "O") {
                    if (verbose_) cout << "previous IC moles for O is: " << icmoles[i] << endl;
                    icmoles[i] = watermoles;
                    solut_->setICmoles(i,(watermoles));
                    if (verbose_) cout << "new ICmoles for O is: " << icmoles[i] << endl;
               }
            }
            double wmv = chemsys_->getNode()->DC_V0(chemsys_->getDCid("H2O@"),
                                                    chemsys_->getP(),
                                                    chemsys_->getTemperature());
            chemsys_->setPhasemass(chemsys_->getPhaseid("aq_gen"),watermass);
            chemsys_->setPhasevolume(chemsys_->getPhaseid("aq_gen"),wmv/watermoles);

            if (verbose_) {
                cout << "Looping over soluble minerals.  ";
                cout << "Here is the list of them:" << endl;
                for (int i = 0; i < solublephase_.size(); i++) {
                    cout << name_[solublephase_[i]] << endl;
                }
                cout.flush();
            }
            for (int i = 0; i < solublephase_.size(); i++) {
                spid = solublephase_[i];
                massdissolved = scaledmass_[spid];
                scaledmass_[spid] = 0.0;
                chemsys_->setMicphasemass(chemsys_->getThermo2mic(spid),
                                      scaledmass_[spid]);
                chemsys_->setMicphasemassdissolved(chemsys_->getThermo2mic(spid),
                                   massdissolved);
                if (verbose_) {
                    cout << name_[spid] << ": Mass dissolved = "
                         << massdissolved << endl;
                    cout.flush();
                    cout << name_[spid] << ": DC molar mass of "
                         << chemsys_->getDCname(chemsysDCid_[spid])
                         << " = " << chemsys_->getDCmolarmass(chemsysDCid_[spid])
                         << endl;
                    cout.flush();
                }
                for (int ii = 0; ii < icmoles.size(); ii++) {
                    icmoles[ii] += ((massdissolved
                           / chemsys_->getDCmolarmass(chemsysDCid_[spid]))
                           * chemsys_->getDCstoich(chemsysDCid_[spid],ii));
                    if (verbose_) {
                        cout << name_[spid] << ":     Total dissolved "
                             << icname[ii] << " = " << icmoles[ii] << " mol "
                             << endl;
                    }
                }
                if (verbose_) cout.flush();

            }

            if (verbose_) {
                cout << "Looping over thermo minerals.  Here is the list of them:" << endl;
                for (int i = 0; i < thermophase_.size(); i++) {
                    cout << name_[thermophase_[i]] << endl;
                }
                cout.flush();
            }
            for (int i = 0; i < thermophase_.size(); i++) {
                tpid = thermophase_[i];
                for (int ii = 0; ii < icmoles.size(); ii++) {
                    icmoles[ii] += ((scaledmass_[tpid]
                                / chemsys_->getDCmolarmass(chemsysDCid_[tpid]))
                                * chemsys_->getDCstoich(chemsysDCid_[tpid],ii));
                    if (verbose_) {
                        cout << name_[tpid] << ":     Total IC " << icname[ii] << " = "
                             << icmoles[ii] << " mol " << endl;
                        cout.flush();
                    }
                }
            }

            // Modify initial pore solution composition if desired
            // input units are mol/kgw
            // watermass is in units of grams, not kg
            
            // icmoles[chemsys_->getICid("K")] += (1.00 * watermass * 0.001);
            // icmoles[chemsys_->getICid("O")] += (1.00 * watermass * 0.001);
            // icmoles[chemsys_->getICid("H")] += (1.00 * watermass * 0.001);

            double kgwatermass = watermass / 1000.0;
            
            map<int,double> iscomp = chemsys_->getInitialSolutionComp();
            map<int,double>::iterator p = iscomp.begin();

            while (p != iscomp.end()) {
                if (verbose_) {
                    cout << "--->Adding " << p->second << " mol/kgw of "
                         << icname[p->first] << " to initial solution." << endl;
                    cout.flush();
                }
                icmoles[p->first] += (p->second * kgwatermass);
                p++;
            }

        }   // End of special first-time tasks

        vector<int> phaselist;
        double icmol,icstoich;
        double phmass,dphmass,partitionedmoles,denom;
        double partitionedmoles_K = 0.0, partitionedmoles_Na = 0.0;
        if (verbose_) cout << "Aqueous mass = " << chemsys_->getPhasemass("aq_gen") << endl;
        for (int i = 0; i < thermophase_.size(); i++) {
            tpid = thermophase_[i];
            if (verbose_) cout << "Partitioning alkali for phase " << name_[tpid] << ": ";
            phaselist = chemsys_->getMic2phase(chemsys_->getThermo2mic(tpid));
            // List of all phases making up this microstructure phase
            // Compute the mass of each phase and sum to get phase mass
            phmass = dphmass = 0.0;
            for (int j = 0; j < phaselist.size(); j++) {
                // Mass will be in g, not kg
                if (verbose_) {
                    cout << endl;
                    cout << "phasename: " << chemsys_->getPhasename(phaselist[j]) << endl;
                    cout << "phasemass: " << chemsys_->getPhasemass(phaselist[j]) << endl;
                }
                dphmass += (double)(chemsys_->getPhasemass(phaselist[j])
                                  - chemsys_->getOphasemass(phaselist[j]));
                phmass += (double)(chemsys_->getPhasemass(phaselist[j]));
            }
            if (verbose_) {
                cout << "Mass is " << phmass << "; mass increment is "
                     << dphmass << endl;
                cout.flush();
            }
            for (int j = 0; j < RdICid_[tpid].size(); j++) {
                icid = RdICid_[tpid][j];
                Rd = Rd_[tpid][j];
                icmol = chemsys_->getVphasestoich("aq_gen",icid);
                if (verbose_) {
                    cout << " Partitioning IC "
                         << chemsys_->getICname(icid) << " with aqueous moles = ";
                    cout << icmol << endl;
                    cout << " Partitioning IC " << chemsys_->getICname(icid);
                    cout << " with actual aqueous moles = ";
                    cout << chemsys_->getPhasestoich(chemsys_->getPhaseid("aq_gen"),icid) << endl;
                    cout << "moles of 'O' in the solution: ";
                    cout << chemsys_->getPhasestoich(chemsys_->getPhaseid("aq_gen"),
                                                     chemsys_->getICid("O")) << endl;
                    cout.flush();
                }
                if (dphmass > 0.0 && Rd > 0.0 && icmol > 0.0) {
                    // phase volume will be in m3, so transform to cm3
                    denom = 1.0 + (1.0e6 * chemsys_->getPhasevolume("aq_gen")
                                           / dphmass / Rd);
                    partitionedmoles = (icmol/denom);
            
                    // be careful whether partitionedmoles is the normalized ones?
                    icmol -= (partitionedmoles);
                    if (icid == chemsys_->getICid("Na")) {
                        partitionedmoles_Na = partitionedmoles * 
                                      chemsys_->getPhasestoich(chemsys_->getPhaseid("aq_gen"),
                                      chemsys_->getICid("O"));
                        if (verbose_) {
                            cout << "partitionedmoles for Na: " << partitionedmoles << endl;
                            cout << "actually removed Na: " << partitionedmoles_Na
                                 << " moles" << endl;
                        }
                        icmoles[icid] = icmoles[icid] - partitionedmoles_Na;
                    }
                    if (icid == chemsys_->getICid("K")) {
                        partitionedmoles_K = partitionedmoles * 
                        chemsys_->getPhasestoich(chemsys_->getPhaseid("aq_gen"),
                                         chemsys_->getICid("O"));
                        if (verbose_) {
                            cout << "partitionedmoles for K: " << partitionedmoles << endl;
                            cout << "actually removed K: " << partitionedmoles_K
                                 << " moles" << endl;
                        }
                        icmoles[icid] = icmoles[icid] - partitionedmoles_K;
                    }
                    if (icmol < 1.0e-16) icmol = 1.0e-16;
                    if (verbose_) {
                        cout << "; new conc = " << icmol / chemsys_->getPhasemass("aq_gen")
                              * 1000.0 << " molal" << endl;
                        cout.flush();
                    }
                }
            }
        }

        if (hyd_time >= leach_time_ && hyd_time < sattack_time_) {

            ///
            /// This is the block for simulating leaching kinetics
            ///

            if (verbose_) cout << " begin to leach at time = " << hyd_time << endl;
            double removed_K = 0.0, removed_Na = 0.0, removed_Ca = 0.0, removed_S = 0.0;
            double removed_O = 0.0, removed_H = 0.0;
            int Kid = chemsys_->getDCid("K+");
            int Naid = chemsys_->getDCid("Na+");
            int Caid = chemsys_->getDCid("Ca+2");
            int Sid = chemsys_->getDCid("SO4-2");
	
            // moles of K+, Na+, Ca+2 in solution before 
            double K1_aq, Na1_aq, Ca1_aq, S1_aq;
            // moles of K+, Na+ Ca+2 in solution after 
            double K2_aq, Na2_aq, Ca2_aq, S2_aq;
            double water_mass = 0.0;
            water_mass = chemsys_->getDCmoles(chemsys_->getDCid("H2O@"))
                            * chemsys_->getDCmolarmass("H2O@");
            //mass will be in g, not kg
            K1_aq = (chemsys_->getNode())->Get_cDC((long int) Kid)
                     * water_mass * 1.0e-3; //mass will be in g, not kg
            Na1_aq = (chemsys_->getNode())->Get_cDC((long int) Naid)
                     * water_mass * 1.0e-3; //mass will be in g, not kg
            Ca1_aq = (chemsys_->getNode())->Get_cDC((long int) Caid)
                     * water_mass * 1.0e-3; //mass will be in g, not kg
            S1_aq = (chemsys_->getNode())->Get_cDC((long int) Sid)
                     * water_mass * 1.0e-3;
            K2_aq = ((chemsys_->getNode())->Get_cDC((long int) Kid) * 0.6 > 0.0005) ? 
        	         (chemsys_->getNode())->Get_cDC((long int) Kid) * 0.6 * water_mass * 1.0e-3 : 0.0005 * water_mass * 1.0e-3;
            Na2_aq = ((chemsys_->getNode())->Get_cDC((long int) Naid) * 0.6 > 0.0005) ?
        	         (chemsys_->getNode())->Get_cDC((long int) Naid) * 0.6 * water_mass * 1.0e-3 : 0.0005 * water_mass * 1.0e-3;
            S2_aq = 0.0005 * water_mass * 1.0e-3;
            // Ca2_aq = 0.0005 * water_mass * 1.0e-3;

            double targetconc = 0.0005;		
            if (verbose_) cout << "set target calcium concentration to be: "
                               << targetconc << endl;
    
            if ((chemsys_->getNode())->Get_cDC((long int) Caid) > (targetconc + 0.0001)) {
                Ca2_aq = 1.0e-7 * water_mass * 1.0e-3;
            } else {
                Ca2_aq = targetconc * water_mass * 1.0e-3;
            }
	
            if (verbose_) cout << "molar mass of water is: "
                               << chemsys_->getDCmolarmass("H2O@") << endl;
	
            removed_K = K1_aq - K2_aq;
            removed_Na = Na1_aq - Na2_aq;
            removed_Ca = Ca1_aq - Ca2_aq;
            if (S1_aq > S2_aq){
                removed_S = S1_aq - S2_aq;
            } else {
                removed_S = 0.0;
            }

            for (int ii=0; ii < icmoles.size(); ii++) {
                if (icname[ii] == "K") {
                    icmoles[ii] = icmoles[ii] - removed_K;
                    if (verbose_) cout << "removed " << removed_K << "moles of K." << endl;
                }
                if (icname[ii] == "Na") {
                    icmoles[ii] = icmoles[ii] - removed_Na;
             	    if (verbose_) cout << "removed " << removed_Na << "moles of Na." << endl;
                }
                if (icname[ii] == "Ca") {
                    icmoles[ii] = icmoles[ii] - removed_Ca;
             	    if (verbose_) cout << "removed " << removed_Ca << "moles of Ca." << endl;
                }
                if (icname[ii] == "S") {
                    icmoles[ii] = icmoles[ii] - removed_S;
                    if (verbose_) cout << "removed " << removed_S << "moles of S." << endl;
                }
                if (icname[ii] == "O") {
                    removed_O = removed_Na + removed_K + 2 * removed_Ca + 4 * removed_S;
                    icmoles[ii] = icmoles[ii] - removed_O;
             	    if (verbose_) cout << "removed " << removed_O << "moles of O." << endl;
                }
                if (icname[ii] == "H") {
                    removed_H = removed_Na + removed_K + 2 * removed_Ca;
                    icmoles[ii] = icmoles[ii] - removed_H;
             	    if (verbose_) cout << "removed " << removed_H << "moles of H." << endl;
                }
            }

        }   // End of block on leaching kinetics

        for (int ii = 0; ii < icmoles.size(); ii++) {
            chemsys_->setICmoles(ii,icmoles[ii]);
            if (verbose_) cout << "ICmoles of " << icname[ii] << " is: " << icmoles[ii] << endl;
        }

        if (hyd_time >= sattack_time_ && hyd_time < leach_time_) {
    
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
 
            int Sid = chemsys_->getDCid("SO4-2");
            int Naid = chemsys_->getDCid("Na+");
            int Kid = chemsys_->getDCid("K+");
            int Caid = chemsys_->getDCid("Ca+2");
            int Caid1 = chemsys_->getDCid("CaOH+");
            int Hid = chemsys_->getDCid("H+");
            int OHid = chemsys_->getDCid("OH-");


            // moles of components in solution before 
            double Na1_aq, K1_aq, Mg1_aq, S1_aq, Ca1_aq, H1_aq, OH1_aq;

            // moles of solution in solution after
            double Na2_aq, K2_aq, Mg2_aq, S2_aq, Ca2_aq, H2_aq, OH2_aq;

            double water_mass = 0.0;
            water_mass = chemsys_->getDCmoles(chemsys_->getDCid("H2O@"))
                       * chemsys_->getDCmolarmass("H2O@");


            double c1; //current concentration
            double c2; //new concentration

            ifstream in("targetConc.dat");
            double concbuff = 0.0;
            in >> concbuff;
            double Na_target = concbuff;
            if (verbose_) cout << "Na_target = " << Na_target << endl;
            in >> concbuff;
            double S_target = concbuff;
            if (verbose_) cout << "S_target = " << S_target << endl;
            in >> concbuff;
            double Ca_target = concbuff;
            if (verbose_) cout << "Ca_target = " << Ca_target << endl;
            in >> concbuff;
            double Ca1_target = concbuff;
            if (verbose_) cout << "Ca1_target = " << Ca1_target << endl;
            in >> concbuff;
            double K_target = concbuff;
            if (verbose_) cout << "K_target = " << K_target << endl;
            in >> concbuff;

            in.close();    

            //double Na_target = 0.15313262, S_target = 0.00196016;
            //double Ca_target = 0.00012733, Ca1_target = 0.00034906;
            //double K_target = 0.4069377;
    
            for (int j = 0; j < 4; j++) { 
              switch(j) {
                case 0: // Na+
                  c1 = (chemsys_->getNode())->Get_cDC((long int) Naid);
                  c2 = Na_target;
                  Na1_aq = c1 * water_mass * 1.0e-3;
                  Na2_aq = c2 * water_mass * 1.0e-3;
                  add_Na += (Na2_aq - Na1_aq);
                  break;

                case 1: // SO4-2
                  c1 = (chemsys_->getNode())->Get_cDC((long int) Sid);
                  c2 = S_target;
                  S1_aq = c1 * water_mass * 1.0e-3;
                  S2_aq = c2 * water_mass * 1.0e-3;
                  add_S = S2_aq - S1_aq;
                  add_O += (add_S * 4.0);
                  break;

                case 2: // K+
                  c1 = (chemsys_->getNode())->Get_cDC((long int) Kid);
                  c2 = K_target;
                  if (c1 >= 0.0005) {
                    K1_aq = c1 * water_mass * 1.0e-3;
                    K2_aq = c2 * water_mass * 1.0e-3;
                    remove_K += (K1_aq - K2_aq);
                  } else {
                    remove_K = 0.0;
                  }
                  break;

                case 3: // Ca2+ and CaOH+  
                  c1 = ((chemsys_->getNode())->Get_cDC((long int) Caid)+ 
                           (chemsys_->getNode())->Get_cDC((long int) Caid1));
                  c2 = Ca_target + Ca1_target;
                  Ca1_aq = c1 * water_mass * 1.0e-3;
                  Ca2_aq = c2 * water_mass * 1.0e-3;
                  remove_Ca += (Ca1_aq - Ca2_aq);
                  remove_O += (remove_K + remove_Ca * 2);
                  remove_H += (remove_K + remove_Ca * 2);
                  break;

                /*
                case 4: // H+
                  c1 = ((chemsys_->getNode())->Get_cDC((long int)Hid));        
                  c2 = H_target;
                  H1_aq = c1 * water_mass * 1.0e-3;
                  H2_aq = c2 * water_mass * 1.0e-3;
                  add_H += H2_aq - H1_aq;
                  break;
                */
                /*
                case 4: // OH-
                  c1 = ((chemsys_->getNode())->Get_cDC((long int)OHid));
                  c2 = OH_target;
                  OH1_aq = c1 * water_mass * 1.0e-3;
                  OH2_aq = c2 * water_mass * 1.0e-3;
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
    
            for (int ii=0; ii < icmoles.size(); ii++) {
                if (icname[ii] == "Na") {
                    icmoles[ii] = icmoles[ii] + add_Na;
                    solut_icmoles[ii] = solut_icmoles[ii] + add_Na;
                    if (verbose_) cout << "added " << add_Na << " moles of Na." << endl;
                }
                if (icname[ii] == "S") {
                    icmoles[ii] = icmoles[ii] + add_S;
                    solut_icmoles[ii] = solut_icmoles[ii] + add_S;
                    if (verbose_) cout << "added " << add_S << " moles of S." << endl;
                }
                if (icname[ii] == "K") {
                    icmoles[ii] = icmoles[ii] - remove_K;
                    solut_icmoles[ii] = solut_icmoles[ii] - remove_K;
                    if (verbose_) cout << "removed " << remove_K << " moles of K." << endl;
                }
                if (icname[ii] == "Ca") {
                    icmoles[ii] = icmoles[ii] - remove_Ca;
                    solut_icmoles[ii] = solut_icmoles[ii] - remove_Ca;
                    if (verbose_) cout << "removed " << remove_Ca << " moles of Ca." << endl;
                }
                if (icname[ii] == "O") {
                    icmoles[ii] = icmoles[ii] + add_O - remove_O;
                    solut_icmoles[ii] = solut_icmoles[ii] + add_O - remove_O;
                    if (verbose_) cout << "added " << (add_O - remove_O) << " moles of O." << endl;
                }
                if (icname[ii] == "H") {
                    icmoles[ii] = icmoles[ii] + add_H - remove_H;
                    solut_icmoles[ii] = solut_icmoles[ii] + add_H - remove_H;
                    if (verbose_) cout << "added " << (add_H - remove_H) << " moles of H." << endl;
                }
            }
   
            for (int ii = 0; ii < icmoles.size(); ii++) {
                chemsys_->setICmoles(ii,icmoles[ii]);
                solut_->setICmoles(ii,solut_icmoles[ii]);
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
 
 
void KineticModel::calculatePhaseChange (const int pid,
                                         double k,
                                         double gamma,
                                         double timestep)
{
    ///
    /// Initialize local variables
    ///

    double DCchange = 0.0;
    vector<int> phaseindex;
    phaseindex.clear();
    phaseindex = chemsys_->getMic2phase(pid);

    ///
    /// Should the next block be only for verbose output
    ///

    if (verbose_) {
        cout << "Micphase " << chemsys_->getMicphasename(pid) << " contains phases: "
             << endl;
        for (int i = 0; i < phaseindex.size(); i++) {
            cout << chemsys_->getPhasename(i) << endl;
        }    
    }

    double phasemoles = 0.0;
    double *phasemolesarray = chemsys_->getPhasemoles();
    for (int i = 0; i < phaseindex.size(); i++) {
      phasemoles += phasemolesarray[i];
    }

    vector<double> phasefrac;
    phasefrac.clear();
    phasefrac.resize(phaseindex.size(), 0.0);

    for (int i = 0; i < phaseindex.size(); i++) {
      double A = 0.0; // A is the surface area
      phasefrac[i] = chemsys_->getPhasemoles(phaseindex[i]) / phasemoles;
      A = lattice_->getSurfaceArea(pid) * phasefrac[i];
      double SI = 0.0; // SI is the saturation index for the phase i
      SI = solut_->getSI(phaseindex[i]);
      
      vector<int> dcindex;
      dcindex.clear();
      dcindex = chemsys_->getPhaseDCmembers(phaseindex[i]);
      if (verbose_) cout << "Phase " << chemsys_->getPhasename(phaseindex[i]) 
                         << " contains DC members: " << endl;
      for (int j = 0; j < dcindex.size(); j++) {
        chemsys_->getDCname(dcindex[j]);
      }
      

      vector<double> dcfrac;
      dcfrac.clear();
      dcfrac.resize(dcindex.size(), 0.0);
      double dcmoles = 0.0;
      for (int j = 0; j < dcindex.size(); j++) {
        dcmoles += chemsys_->getDCmoles(dcindex[j]);
      }
      
      double surfaceArea = 0.0;
      for (int j = 0; j < dcindex.size(); j++) {
        dcfrac[j] = chemsys_->getDCmoles(dcindex[j]) / dcmoles;
        surfaceArea = A * dcfrac[j];
        DCchange = k * surfaceArea * pow((SI - 1),gamma);

        ///
        /// Convert time step from days to seconds, because DCchange
        /// has units of mol/s
        ///

        DCchange = DCchange * timestep * 24.0 * 60.0 * 60.0; // DCchange unit: moles/s

        double newDCmoles = chemsys_->getDCmoles(dcindex[j]) + DCchange;
        chemsys_->setDCupperlimit(dcindex[j],newDCmoles);
        chemsys_->setDClowerlimit(dcindex[j],newDCmoles);

      }
    }

    return;
}

void KineticModel::initializeMoles () 
{
    int icnum = chemsys_->getICnum();
    vector<double> icmoles;
    icmoles.clear();
    icmoles.resize(icnum,1.0e-9);
    vector<string> icname = chemsys_->getICname();
    
    try {
        for (int i = 0; i < icmoles.size(); i++) {
            if ((chemsys_->getICname(i)) == "Nit")
                  icmoles[i] = chemsys_->getICmoles(i);
        }
        double watermass = 100.0 * wcratio_;
        int waterid = chemsys_->getDCid("H2O@");
        double watermolarmass = chemsys_->getDCmolarmass(waterid);
        if (watermolarmass <= 0.0) {
            throw FloatException("KineticModel","initializeMoles",
                                 "Divide by zero error");
        }
        for (int ii = 0; ii < icnum; ii++) {
            icmoles[ii] += ((watermass/chemsys_->getDCmolarmass(waterid)
                           * chemsys_->getDCstoich(waterid,ii)));
            if (verbose_) {
                cout << "H2O@ :     Total contributed " << icname[ii] << " = "
                     << icmoles[ii] << " mol " << endl;
                cout.flush();
            }
        }
        for (int ii = 0; ii < icnum; ii++) {
            chemsys_->setICmoles(ii,icmoles[ii]);
        }

        int dcnum = chemsys_->getDCnum();
        vector<double> dcmoles;
        dcmoles.clear();
        dcmoles.resize(dcnum,1.0e-9);
        vector<string> dcname = chemsys_->getDCname();
        for (int i = 0; i < dcmoles.size(); i++) {
            if (dcname[i] == "O2" || dcname[i] == "N2") {
                dcmoles[i] = chemsys_->getDCmoles(i);
            }
        }
        for (int i = 0; i < dcmoles.size(); i++) {
            if (i != waterid) chemsys_->setDCmoles(i,dcmoles[i]);
        }
        chemsys_->setDCmoles(waterid,watermass
                             / chemsys_->getDCmolarmass(waterid));

        if (verbose_) {
            cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
            cout << "REPORTING DC MOLES:" << endl;
            for (int i = 0; i < dcmoles.size(); i++) {
                cout << dcname[i] << ": " << chemsys_->getDCmoles(i) << " mol" << endl;
            }
            cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
            cout.flush();
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
        EOBException ex("KineticModel","setKineticDCmmoles",
                            oor.what(),0,0);
        ex.printException();
        exit(1);
    }
    
}

void KineticModel::setKineticDCmoles ()
{
    if (verbose_) {
        cout << "    In setKineticDCmoles..." << endl;
        cout.flush();
    }

    int kpid,csid;

    try {
        int waterid = chemsys_->getDCid("H2O@");
        double watermoles = chemsys_->getDCmoles(waterid);
        double watermm = chemsys_->getDCmolarmass(waterid);
        double watermass = watermoles * watermm;
        if (verbose_) {
            cout << "        "
                 << chemsys_->getDCname(waterid) << ": Mass = " << watermass
                 << ", Molar mass = " << watermm
                 << endl;
        }
        for (int i = 0; i < kineticphase_.size(); i++) {
            kpid = kineticphase_[i];
            csid = chemsysDCid_[kpid];
            if (chemsys_->getDCmolarmass(csid) <= 0.0) {
                throw FloatException("KineticModel","setKineticDCmoles",
                                     "Divide by zero error");
            }
            if (verbose_) {
                cout << "        Clinker phase "
                     << name_[kpid] << ": Mass = " << scaledmass_[kpid]
                     << ", Molar mass = " << chemsys_->getDCmolarmass(csid)
                     << endl;
            }
            chemsys_->setDCmoles(csid,(scaledmass_[kpid]
                                 / chemsys_->getDCmolarmass(csid)));
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
        EOBException ex("KineticModel","setKineticDCmmoles",
                           oor.what(),0,0);
        ex.printException();
        exit(1);
    }
    return;
}

void KineticModel::zeroKineticDCmoles ()
{
    try {
        for (int i = 0; i < kineticphase_.size(); i++) {
            chemsys_->setDCmoles(chemsysDCid_[kineticphase_[i]],0.0);
        }
    }
    catch (EOBException eex) {
        eex.printException();
        exit(0);
    }
    return;
}
