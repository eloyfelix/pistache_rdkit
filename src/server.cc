#include <pistache/http.h>
#include <pistache/router.h>
#include <pistache/endpoint.h>

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

#include <GraphMol/inchi.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

#include <nlohmann/json.hpp>

using namespace Pistache;
using json = nlohmann::json;


namespace Generic
{
    void handleReady(const Rest::Request &, Http::ResponseWriter response)
    {
        response.send(Http::Code::Ok, "1");
    }

    std::unique_ptr<RDKit::ROMol> readMol(std::string input)
    {
        std::unique_ptr<RDKit::ROMol> mol;
        try
        {
            if (input.find("M  END") != std::string::npos) {
                mol.reset(RDKit::MolBlockToMol(input));
            }
            else
            {
                mol.reset(RDKit::SmilesToMol(input));
            }            
        }
        catch (RDKit::MolSanitizeException &e)
        {
            std::cerr << e.what() << std::endl;
        }
        return mol;
    }
} // namespace Generic

class PistacheRDKit
{
public:
    PistacheRDKit(Address addr) : httpEndpoint(std::make_shared<Http::Endpoint>(addr))
    {
        // setup the PAINS filters
        fcparams.addCatalog(RDKit::FilterCatalogParams::PAINS_A);
        fcparams.addCatalog(RDKit::FilterCatalogParams::PAINS_B);
        fcparams.addCatalog(RDKit::FilterCatalogParams::PAINS_C);
        filterCatalog = fcparams;
    }

    void init(size_t thr = 2)
    {
        auto opts = Http::Endpoint::options().threads(thr);
        httpEndpoint->init(opts);
        setupRoutes();
    }

    void start()
    {
        httpEndpoint->setHandler(router.handler());
        httpEndpoint->serve();
    }

private:
    void setupRoutes()
    {
        using namespace Rest;
        Routes::Get(router, "/ready", Routes::bind(&Generic::handleReady));
        Routes::Post(router, "/mcs", Routes::bind(&PistacheRDKit::mcs, this));
        Routes::Post(router, "/painsFilters", Routes::bind(&PistacheRDKit::painsFilters, this));
        Routes::Post(router, "/molblock2inchi", Routes::bind(&PistacheRDKit::molblock2inchi, this));
        Routes::Post(router, "/mol2inchi", Routes::bind(&PistacheRDKit::mol2inchi, this));
        Routes::Post(router, "/inchi2inchikey", Routes::bind(&PistacheRDKit::inchi2inchikey, this));
        Routes::Post(router, "/descriptors", Routes::bind(&PistacheRDKit::descriptors, this));
    }

    void molblock2inchi(const Rest::Request &request, Http::ResponseWriter response)
    {
        /**
         * Get the InChI for a molblock bypassing RDKit parsing.
         */
        RDKit::ExtraInchiReturnValues tmp;
        std::string inchi = RDKit::MolBlockToInchi(request.body(), tmp);
        response.send(Http::Code::Ok, inchi);
    }

    void mol2inchi(const Rest::Request &request, Http::ResponseWriter response)
    {
        /**
         * Get the InChI for a molblock or SMILES with RDKit parsing.
         */
        std::string input(request.body());
        std::unique_ptr<RDKit::ROMol> mol(Generic::readMol(input));
        if (!mol)
        {
            response.send(Http::Code::Internal_Server_Error, "Cannot create molecule from input");
        }
        else
        {
            RDKit::ExtraInchiReturnValues tmp;
            std::string inchi = RDKit::MolToInchi(*mol, tmp);
            response.send(Http::Code::Ok, inchi);
        }
    }

    void inchi2inchikey(const Rest::Request &request, Http::ResponseWriter response)
    {
        /**
         * Get the InChIKey for an InChI.
         */
        std::string inchi(request.body());
        std::string inchikey = RDKit::InchiToInchiKey(inchi);
        response.send(Http::Code::Ok, inchikey);
    }

    void mcs(const Rest::Request &request, Http::ResponseWriter response)
    {
        /**
         * Find the Maximum Common Substructure of a set of SMILES.
         */
        std::string sbody(request.body());

        // split SMILES in body
        std::stringstream ss(sbody);
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> smiles(begin, end);
        std::copy(smiles.begin(), smiles.end(), std::ostream_iterator<std::string>(std::cout, "\n"));

        std::vector<RDKit::ROMOL_SPTR> mols;
        for(const auto& sm: smiles) {

            RDKit::ROMOL_SPTR mol;
            try
            {
                mol.reset(RDKit::SmilesToMol(sm));
                if (mol){
                    mols.emplace_back(mol);
                }
                else
                {
                    std::cerr << "Cannot create molecule from : '" << sm << "'" << std::endl;
                }
            }
            catch (RDKit::MolSanitizeException &e)
            {
                std::cerr << e.what() << "\t" << sm << std::endl;
            }
        }
        RDKit::MCSResult res = RDKit::findMCS(mols);
        response.send(Http::Code::Ok, res.SmartsString);
    }

    void painsFilters(const Rest::Request &request, Http::ResponseWriter response)
    {
        /**
         * Get PAINS filters alerts for a compound.
         */
        std::string input(request.body());
        std::unique_ptr<RDKit::ROMol> mol(Generic::readMol(input));

        std::vector<std::string> alerts;
        if (!mol)
        {
            response.send(Http::Code::Internal_Server_Error, "Cannot create molecule from input");
        }
        else
        {
            const RDKit::FilterCatalog::CONST_SENTRY entry = filterCatalog.getFirstMatch(*mol);
            if (entry) {
                alerts.push_back(entry->getDescription());
            }   
        }
        json out(alerts);
        response.send(Http::Code::Ok, out.dump());
    }

    void descriptors(const Rest::Request &request, Http::ResponseWriter response)
    {
        /**
         * Get a set of common descriptors for a compound.
         */
        std::string input(request.body());
        std::unique_ptr<RDKit::ROMol> mol(Generic::readMol(input));
        std::map<std::string, double> res_map;
        if (!mol)
        {
            response.send(Http::Code::Internal_Server_Error, "Cannot create molecule from input");
        }
        else
        {
            res_map["logP"] = RDKit::Descriptors::calcClogP(*mol);;
            res_map["mw"] = RDKit::Descriptors::calcExactMW(*mol);
            res_map["rtb"] = RDKit::Descriptors::calcNumRotatableBonds(*mol);
            res_map["hba"] = RDKit::Descriptors::calcNumHBA(*mol);
            res_map["hbd"] = RDKit::Descriptors::calcNumHBD(*mol);
            res_map["tpsa"] = RDKit::Descriptors::calcTPSA(*mol);
            res_map["nrings"] = RDKit::Descriptors::calcNumRings(*mol);
            res_map["nha"] = mol->getNumHeavyAtoms();

        }
        json out(res_map);
        response.send(Http::Code::Ok, out.dump());
    }


    Rest::Router router;
    std::shared_ptr<Http::Endpoint> httpEndpoint;
    RDKit::FilterCatalogParams fcparams;
    RDKit::FilterCatalog filterCatalog;
};

int main(int argc, char *argv[])
{
    Port port(9080);

    int thr = 2;

    if (argc >= 2)
    {
        port = std::stol(argv[1]);
        if (argc == 3)
            thr = std::stol(argv[2]);
    }

    Address addr(Ipv4::any(), port);

    PistacheRDKit pr(addr);

    pr.init(thr);
    std::cout << "Pistache RDKit API started" << std::endl;
    pr.start();
}
