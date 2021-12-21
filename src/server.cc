#include <nlohmann/json.hpp>

#include <pistache/http.h>
#include <pistache/router.h>
#include <pistache/endpoint.h>

#include <GraphMol/inchi.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/MolHash/MolHash.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

using namespace Pistache;
using json = nlohmann::json;

namespace Generic
{
    void handleReady(const Rest::Request &, Http::ResponseWriter response)
    {
        response.send(Http::Code::Ok, "1");
    }

    template<typename MOLTYPE>
    std::unique_ptr<MOLTYPE> readMol(std::string input)
    {
        std::unique_ptr<MOLTYPE> mol;
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
        Routes::Post(router, "/murckoScaffold", Routes::bind(&PistacheRDKit::murckoScaffold, this));
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
        std::unique_ptr<RDKit::ROMol> mol(Generic::readMol<RDKit::ROMol>(input));
        if (!mol)
        {
            response.send(Http::Code::Internal_Server_Error, "Can't create mol object from input");
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
        std::stringstream ss(sbody);
        std::vector<std::string> smiles;
        std::string token;
        while (std::getline(ss, token, '\n')) {
            smiles.push_back(token);
        }

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
                    std::cerr << "Can't create mol object from : '" << sm << "'" << std::endl;
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
        std::unique_ptr<RDKit::ROMol> mol(Generic::readMol<RDKit::ROMol>(input));
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
        std::unique_ptr<RDKit::ROMol> mol(Generic::readMol<RDKit::ROMol>(input));
        std::map<std::string, double> res_map;
        if (!mol)
        {
            response.send(Http::Code::Internal_Server_Error, "Can't create mol object from input");
        }
        else
        {
            res_map["ClogP"] = RDKit::Descriptors::calcClogP(*mol);
            res_map["ExactMW"] = RDKit::Descriptors::calcExactMW(*mol);
            res_map["NumRotatableBonds"] = RDKit::Descriptors::calcNumRotatableBonds(*mol);
            res_map["NumHBA"] = RDKit::Descriptors::calcNumHBA(*mol);
            res_map["NumHBD"] = RDKit::Descriptors::calcNumHBD(*mol);
            res_map["TPSA"] = RDKit::Descriptors::calcTPSA(*mol);
            res_map["NumRings"] = RDKit::Descriptors::calcNumRings(*mol);
            res_map["NumHeavyAtoms"] = mol->getNumHeavyAtoms();
        }
        json out(res_map);
        response.send(Http::Code::Ok, out.dump());
    }

    void murckoScaffold(const Rest::Request &request, Http::ResponseWriter response)
    {
        /**
         * Get the Bemis-Murcko scaffold for a molecule.
         */
        std::string input(request.body());
        std::unique_ptr<RDKit::RWMol> mol(Generic::readMol<RDKit::RWMol>(input));
        std::string scaffold;
        if (!mol)
        {
            response.send(Http::Code::Internal_Server_Error, "Can't create mol object from input");
        }
        else
        {
            scaffold = RDKit::MolHash::MolHash(mol.get(), RDKit::MolHash::HashFunction::MurckoScaffold);
        }
        response.send(Http::Code::Ok, scaffold);
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
