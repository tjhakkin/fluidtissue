//
// Returns a subset of given objects from a given Wavefront OBJ file.
// Writes output file as [input name]_[object 1]_[object 2].obj
//
// Usage:
// parseobj [input file] [object 1] [object 2] ...
//

#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include "tinyobjloader/tiny_obj_loader.h"



int parseObjData( const std::vector<std::string>& args, std::ofstream& out )
{
    tinyobj::attrib_t                   attrib;
    std::vector<tinyobj::shape_t>       shapes;
    std::vector<tinyobj::material_t>    materials;
    std::string                         error;
    bool ret = tinyobj::LoadObj( &attrib, &shapes, &materials, &error, args[1].c_str() );
    if (!ret) {
        std::cerr << error << std::endl;
        return -1;
    }

    std::cout << "# of vertices  = " << attrib.vertices.size() / 3 << std::endl;
    std::cout << "# of normals   = " << attrib.normals.size() / 3 << std::endl;
    std::cout << "# of texcoords = " << attrib.texcoords.size() / 2 << std::endl;
    std::cout << "# of materials = " << materials.size() << std::endl;
    std::cout << "# of shapes    = " << shapes.size() << std::endl;
    std::cout << "Found shapes:";
    for (size_t i=0; i<shapes.size(); i++)
        std::cout << " " << shapes[i].name;
    std::cout << std::endl << std::endl;

    std::vector<std::string> outputShapes;          // shapes to be extracted
    std::copy(args.begin()+2, args.end(), std::back_inserter(outputShapes));
    std::cout << "Extracting shapes:"; 
    for (size_t i=0; i<outputShapes.size(); i++)
        std::cout << " " << outputShapes[i];
    std::cout << std::endl;

    // face indexing offset; Wavefront indexing begins at 1.
    size_t faceOffset = 1;

    for (size_t i=0; i<shapes.size(); i++) {    // loop over objects
        if (std::find(outputShapes.begin(), outputShapes.end(), shapes[i].name) == outputShapes.end())
            continue;

        size_t nt = shapes[i].mesh.num_face_vertices.size();    // # of faces
        std::vector<int32_t> vIdx(3*nt);       // vertex indices
        std::vector<int32_t> nIdx(3*nt);       // vertex normal indices

        // read vertex, normal global indices.
        for (size_t j=0; j<nt; j++) {          // loop over object faces
            int n = shapes[i].mesh.num_face_vertices[j];
            if (n != 3) {
                std::cerr << __FUNCTION__ << ": encoutered face with " << n
                          << " vertices. Only triangles supported!." << std::endl;
                return -1;
            }

            for (int k=0; k<3; k++) {       // loop over triangle vertices
                vIdx[3*j+k] = shapes[i].mesh.indices[3*j+k].vertex_index;
                nIdx[3*j+k] = shapes[i].mesh.indices[3*j+k].normal_index;
            }
        }

        // construct map from scene-global vertex indices to object-local indices
        std::vector<int32_t> vIdxSort = vIdx;
        std::sort( vIdxSort.begin(), vIdxSort.end() );
        auto last = std::unique( vIdxSort.begin(), vIdxSort.end() );
        vIdxSort.erase( last, vIdxSort.end() );
        size_t np = vIdxSort.size();
        std::map<int,int> idxMap;
        for (size_t j=0; j<np; j++)
            idxMap[vIdxSort[j]] = j;

        // write object vertices
        out << std::endl;
        for (size_t j=0; j<np; j++) {
            uint32_t k = 3 * vIdxSort[j];      // vertex position start index
            out << "v " << attrib.vertices[k+0] << " " 
                        << attrib.vertices[k+1] << " "
                        << attrib.vertices[k+2] << std::endl;
        }
        out << std::endl;

        // write vertex normals
        for (size_t j=0; j<np; j++) {
            int32_t k = 3 * nIdx[j];      // vertex position start index
            if (k < 0)
                continue;
            out << "vn " << attrib.normals[3*j+0] << " " 
                         << attrib.normals[3*j+1] << " "
                         << attrib.normals[3*j+2] << std::endl;
        }
        out << std::endl;

        out << "g " << shapes[i].name << std::endl;
        // NOTE: Using the material of the first face for the whole shape!
        int matId = shapes[i].mesh.material_ids[0];
        out << "usemtl " << materials[matId].name << std::endl;

        // write faces
        for (size_t j=0; j<nt; j++) {
            // Wavefront indexing begins at 1, hence +1.
            out << "f " << idxMap[vIdx[3*j+0]] + faceOffset << "//";
            if (nIdx[3*j+0] >= 0)
                out << idxMap[nIdx[3*j+0]] + faceOffset;
            out << " " << idxMap[vIdx[3*j+1]] + faceOffset << "//";
            if (nIdx[3*j+1] >= 0)
                out << idxMap[nIdx[3*j+1]] + faceOffset;
            out << " " << idxMap[vIdx[3*j+2]] + faceOffset << "//";
            if (nIdx[3*j+2] >= 0)
                out << idxMap[nIdx[3*j+2]] + faceOffset;
            out << std::endl;
        }

        faceOffset += np;
    }

    return 0;
}



int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv+argc);
    if (argc < 3) {
        std::cout << "Usage: parseobj [input file] [object 1] [object 2] ..." << std::endl;
        return 0;
    }

    // Form output file name as [input]_[output shapes].obj
    std::string input = args[1];
    size_t idx = input.find_last_of(".");
    if (idx == std::string::npos)
        return -1;   
    std::string output = input.substr(0, idx);
    for (size_t i=2; i<args.size(); i++)
        output = output + "_" + args[i];
    output = output + ".obj";

    // Open output file
    std::ofstream out(output);
    if (!out.good()) {
        std::cerr << "Error: Cannot open " << output << " for writing."
                  << std::endl;
        return -1;
    }

    // find materials file definition, write to output.
    std::ifstream in(input);
    std::string s;
    std::string mtllib = "";
    in >> s;
    if (!s.find("mtllib"))
        in >> mtllib;
    out << "mtllib " << mtllib << std::endl;

    // Parse & write geometry data.
    parseObjData(args, out);

    out.close();

    return 0;
}
