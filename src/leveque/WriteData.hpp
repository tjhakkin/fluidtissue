//
// Writes model geometry output to Wavefront OBJ file, additional data to .txt files.
//

#pragma once

#include "fem/MathTypes.hpp"
#include "fem/LevelSet.hpp"


namespace cfish {

// Types of data to be exported by Write_*() functions.
enum { INTERIOR = 0x01, INTERFACE = 0x02, BOUNDARY = 0x04, KNOTS = 0x08 };


//
// Writes the requested objects to a Wavefront OBJ file.
//
int     Write_Wavefront_obj( const LevelSet& lset, const std::string& output_fname, 
                             const int exportData );
}
