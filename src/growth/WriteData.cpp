#include <fstream>

#include "WriteData.hpp"
#include "common/Utils.hpp"

using namespace FEM;


int cfish::Write_Wavefront_obj( const LevelSet& lset, const RD& rd, const Knots& knots,
                                const std::string& output_fname, const int exportData )
{
    //
    // Write .obj file (geometry)
    //
    std::ofstream out( output_fname + ".obj" );
    if (!out.good())
        return -1;

    // Set materials reference; the materials file is written after the OBJ.
    out << "mtllib " << output_fname << ".mtl" << "\n\n";

    // for OBJ indices must start from 1
    size_t offset = 1;
    auto& mesh = lset.getMesh();

    if (exportData & INTERIOR) {
        out << "#\n";
        out << "# shape interior\n";
        out << "#\n\n";

        auto& meshRD = rd.getMesh();
        for (auto p : meshRD.p())
            out << "v " << p(0) << " " << p(1) << " " << p(2) << "\n";
        out << "\n";

        out << "g interior\n";
        for (Vec4i t : meshRD.t()) {
            t += offset * Vec4i::Ones();
            // tetrahedra to triangles:
            out << "f " << t(2) << " " << t(1) << " " << t(0) << "\n";
            out << "f " << t(1) << " " << t(2) << " " << t(3) << "\n";
            out << "f " << t(0) << " " << t(3) << " " << t(2) << "\n";
            out << "f " << t(3) << " " << t(0) << " " << t(1) << "\n";               
        }
        offset += meshRD.p().size();
    }

    if (exportData & INTERFACE) {
        out << "#\n";
        out << "# shape interface\n";
        out << "#\n\n";

        // Interface mesh (implicit surface reconstruction).
        auto& ifMesh = lset.getInterface();
        for (auto& p : ifMesh.p())                        // vertices
            out << "v " << p(0) << " " << p(1) << " " << p(2) << "\n";
        out << "\n";

        // Compute & write vertex normals
        for (size_t i=0; i<ifMesh.p().size(); i++) {
            Vec3f normal(0.0f, 0.0f, 0.0f);
            auto& p2t = ifMesh.p2t(i);
            for (size_t j=0; j<p2t.size(); j++) {
                auto& t = ifMesh.t(p2t(j));
                Vec3f e1 = ifMesh.p(t(1)) - ifMesh.p(t(0));
                Vec3f e2 = ifMesh.p(t(2)) - ifMesh.p(t(0));
                normal += e1.cross(e2);
            }
            out << "vn " << normal(0) << " " << normal(1) << " " << normal(2) << "\n";
        }
        out << "\n";

        out << "g interface\n";
        out << "usemtl 0\n";
        for (Vec4i t : ifMesh.t()) {                      // triangles
            t += offset * Vec4i::Ones();
            out << "f " << t(0) << "//" << t(0) << " "  // with normals
                        << t(1) << "//" << t(1) << " "
                        << t(2) << "//" << t(2) << "\n";
        }
        offset += ifMesh.p().size();
    }

    if (exportData & BOUNDARY) {
        out << "#\n";
        out << "# domain boundary\n";
        out << "#\n\n";

        // Get unique boundary nodes
        std::vector<int> nodes = unique_nodes( mesh.b_tris() );
        // Node map from global node index to boundary-local index.
        std::map<int,int> imap;
        for (size_t i=0; i<nodes.size(); i++)
            imap[ nodes[i] ] = i;

        for (auto& node : nodes) {
            auto& p = mesh.p( node );
            out << "v " << p(0) << " " << p(1) << " " << p(2) << "\n";
        }
        out << "\n";

        out << "g boundary\n";
        out << "usemtl 1\n";
        for (auto& t : mesh.b_tris()) {
            out << "f " << offset + imap[t(0)] << " " << offset + imap[t(1)] << " "
                        << offset + imap[t(2)] << "\n";
        }
        offset += nodes.size();
    }

    if (exportData & KNOTS) {
        out << "#\n";
        out << "# Center knots\n";
        out << "#\n\n";
        std::vector<Vec3f> vertices;
        std::vector<Vec3i> triangles;
        knots.knotsToMesh( vertices, triangles );
        for (auto& p : vertices)
            out << "v " << p(0) << " " << p(1) << " " << p(2) << "\n";
        out << "\n";

        out << "g knots\n";
        out << "usemtl 2\n";
        for (Vec3i t : triangles) {                     // triangles
            t += offset * Vec3i::Ones();
            out << "f " << t(0) << " " << t(1) << " " << t(2) << "\n";
        }
    }
    out.close();

    //
    // Write .mtl file (materials)
    //
    out.open( output_fname + ".mtl" );
    if (!out.good())
        return -1;

    // Diffuse interface in gray
    out << "newmtl 0" << "\n";
    out << "Ka 0 0 0" << "\n";
    out << "Kd 0.6 0.6 0.6" << "\n";
    out << "d 1" << "\n";
    out << "Ks 0 0 0" << "\n";
    out << "Ns 0 interface" << "\n\n";

    // Transparent boundary in blue
    out << "newmtl 1" << "\n";
    out << "Ka 0 0 0" << "\n";
    out << "Kd 0.0 0.0 1.0" << "\n";
    out << "d 0.2" << "\n";
    out << "Ks 0 0 0" << "\n";
    out << "Ns 0 boundary" << "\n";

    // Knots in red
    out << "newmtl 2" << "\n";
    out << "Ka 0 0 0" << "\n";
    out << "Kd 1.0 0.0 0.0" << "\n";
    out << "d 1" << "\n";
    out << "Ks 0 0 0" << "\n";
    out << "Ns 0 knots" << "\n";

    out.close();

    return 0;
}



int cfish::Write_data( const LevelSet& lset, const RD& rd, const VecXd& flow, 
                       const std::string& output_fname, 
                       const int exportData )
{
    auto& mesh = lset.getMesh();
    size_t n = mesh.p().size();
    std::ofstream out;

    //
    // Data for interior element visualization.
    //
    if (exportData & INTERIOR) {
        out.open( output_fname + "_interior.txt" );
        if (!out.good())
            return -1;

        const VecXd& act = rd.getActivator();
        const VecXd& gf = rd.getGrowthFactor();
        const std::vector<int>& rd_nodes = rd.getNodes();

        out << "Activator\tGrowth factor\n";
        for (int i : rd_nodes)
            out << act(i) << "\t" << gf(i) << "\t" << "\n";
        out.close();
    }

    //
    // Full domain pressure with mesh coordinates (not used for visualization).
    //
    out.open( output_fname + "_pressure.txt" );
    if (!out.good())
        return -1;

    out << "X\tY\tZ\tPressure\n";
    VecXd p = flow.tail(n);     // pressure
    for (size_t i=0; i<n; i++)
        out << mesh.p(i)(0) << "\t" << mesh.p(i)(1) << "\t" << mesh.p(i)(2) 
            << "\t" << p(i) << "\n";
    out.close();

    return 0;
}
