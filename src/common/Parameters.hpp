//
// A representation of model parameters. Initialize and read a parameters in
// file "par.txt" in ToothMaker/mbaker format with
//
// Parameters par;
// par.import("par.txt");
//
// Access the value (type double) of parameter "Act" with
//
// par["Act"];
//

#pragma once

#include <iostream>
#include <fstream>
#include <map>


struct Parameters {
    typedef std::map<std::string, double> parameterMap;


    double operator[]( const std::string& key ) const
    { 
        try {
            return m_variables.at(key);
        }
        catch (const std::exception& e) {
            std::cerr << "** Warning: Parameter '" << key << "' not defined. "
                      << "Defaulting to value 0.0." << std::endl; 
            return 0.0;
        }
    } 


    int import( const std::string& pfile )
    {
        std::ifstream in( pfile );

        while (true) {
            if (!in.good())
                return -1;
            in >> std::ws;           
            if (in.eof())
                return 0;
            
            std::string line;
            std::getline( in, line );
            if (line[0] == '#')
                continue;
            std::string name = line.substr( 0, line.find("==") );
            if (name == "model" || name == "viewthresh" ||
                name == "viewmode" || name == "iter")
                continue;           // ignore mbaker internal variable names
            std::string value = line.substr( line.find("==")+2 );

            try {
                m_variables[ name ] = std::stod( value );
            }
            catch (std::exception& e) {
                std::cerr << "Exception: " << e.what() << " for value '" << value 
                        << "'." << std::endl;
                return -1;
            }
        }

        return 0;
    }


    const parameterMap& variables() const       { return m_variables; };


    private:
        parameterMap  m_variables;
};
