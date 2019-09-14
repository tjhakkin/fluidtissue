//
// Shared main() and argument handler for both growth and leveque binaries; 
// select by setting -DCFISH_LEVEQUE or -DCFISH_GROWTH compiler flag.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

#ifdef CFISH_LEVEQUE
#include "leveque/leveque.hpp"
#endif
#ifdef CFISH_GROWTH
#include "growth/growth.hpp"
#endif

#include "common/Parameters.hpp"
#include "common/Utils.hpp"



void printHelp( const std::string& name )
{
    std::cout << "Usage: " << name << " --param [parameters] --niter [iterations] "
                 "--step [output stepsize] --id [run id]" << std::endl;
    return;
}



int parseArguments( const std::vector<std::string>& args, std::string& pfile,
                    int& iterations, int& step, int&id )
{
    if (args.size() < 8)
        return 0;

    int argsFound = 0;

    auto res = std::find( std::begin(args), std::end(args), "--param");
    if (res != std::end(args)) {
        pfile = *(res+1);
        argsFound++;    
    }
    res = std::find( std::begin(args), std::end(args), "--niter");
    if (res != std::end(args)) {
        iterations = std::stoi( *(res+1) );
        argsFound++;
    }
    res = std::find( std::begin(args), std::end(args), "--step");
    if (res != std::end(args)) {
        step = std::stoi( *(res+1) );
        argsFound++;
    }
    res = std::find( std::begin(args), std::end(args), "--id");
    if (res != std::end(args)) {
        id = std::stoi( *(res+1) );
        argsFound++;
    }

    return argsFound;
}



int main( int argc, char** argv )
{
    FEM::chrono_time  time_start;
    FEM::tic( time_start );

    std::string pfile;          // parameters file
    int iterations, step, id;   // nof. iterations, results step size, run id
    std::vector<std::string> args( argv+1, argv+argc );
    int n = parseArguments( args, pfile, iterations, step, id );
    if (n < 4 || args[0] == "--help") {
        printHelp( argv[0] );
        return 0;
    }

    Parameters parameters;
    if (parameters.import(pfile)) {
        std::cerr << "Error: Failed to import parameters from '" << pfile << "'." 
                  << std::endl;
        return -1;
    }

    std::cout << "** Parameters:" << std::endl;
    for (auto p : parameters.variables()) {
         std::cout << std::setw(20) << std::left << p.first << ": " << p.second 
                   << std::endl;
    }

    #ifdef CFISH_LEVEQUE
    LeVeque app;
    #endif
    #ifdef CFISH_GROWTH
    Growth app;
    #endif

    if (app.init( parameters, id ))
        return -1;
    app.run( iterations );

    std::cout << "Total execution time: " << FEM::toc( time_start ) / 1000.0f 
              << " s." << std::endl;

    return 0;
}
