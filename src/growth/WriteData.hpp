//
// Writes model geometry output to Wavefront OBJ file, additional data to .txt files.
//

#pragma once

#include "Knots.hpp"
#include "RD.hpp"
#include "fem/MathTypes.hpp"
#include "fem/LevelSet.hpp"


namespace cfish {

// Types of data to be exported by Write_*() functions.
enum { INTERIOR = 0x01, INTERFACE = 0x02, BOUNDARY = 0x04, KNOTS = 0x08 };


//
// Writes the requested objects to a Wavefront OBJ file.
//
int     Write_Wavefront_obj( const LevelSet& lset, const RD& rd, const Knots& knots,
                             const std::string& output_fname, const int exportData );

//
// Writes morphogen concentrations, Stokes pressure into tab-separated txt files.
//
int     Write_data( const LevelSet& lset, const RD& rd, const FEM::VecXd& flow, 
                    const std::string& output_fname, const int exportData );
   
}
