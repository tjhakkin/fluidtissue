#include <iostream>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#include "growth.hpp"
#include "WriteData.hpp"
#include "fem/TetraMesh.hpp"
#include "common/Utils.hpp"

using namespace FEM;


int Growth::init( const Parameters& parameters, const int id )
{
    m_parameters = parameters;
    m_id = id;

    chrono_time timeInit; 
    std::cout << std::endl << "** Initializing..." << std::endl; tic( timeInit );

    // Assuming hyperthreading is enabled, often the best performance with Eigen
    // is obtained by using the physical number of cores, or nof. threads / 2.
    int nThreads = std::thread::hardware_concurrency() / 2;
    std::cout << "Number of threads: " << nThreads << std::endl;
    Eigen::setNbThreads( nThreads );


    chrono_time timeStart;
    std::cout << "1/3 : Initializing level set..." << std::endl; tic( timeStart );
    if (m_lset.init( m_parameters ))
        return -1;
    std::cout << "  Mesh size: " << m_lset.getMesh().sizeInMemory() << " MB." << std::endl;
    std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;


    std::cout << "2/3 : Initializing solvers..." << std::endl; tic( timeStart );
    size_t n = m_lset.getMesh().p().size();
    m_flow.setZero(4*n);        // 3D velocity + pressure
    m_rd.init( m_lset, m_parameters );
    m_knots.init( m_parameters );
    if (m_stokes.init( m_lset, m_parameters )) {
        std::cerr << "Error: Stokes velocity degree must be >= 1." << std::endl;
        return -1;
    }
    std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;


    std::cout << "3/3 : Storing initial state." << std::endl; tic( timeStart );
    if (storeModelState_(0) < 0) {
        std::cerr << "Error: Cannot write model state to output files." << std::endl;
        return -1;
    }
    std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;

    std::cout << "Initialization took " << toc( timeInit ) << " ms." << std::endl;

    return 0;
}



int Growth::run( const int iterations )
{
    for (int i=1; i<iterations+1; i++) {
        chrono_time timeIter;
        std::cout << std::endl << "Iteration " << i << std::endl; tic( timeIter );


        chrono_time timeStart;
        std::cout << "* Solving Stokes." << std::endl; tic( timeStart );
        // By default, growth factor is the source morphogen for if the Stokes.
        // Set activator as the source for free-form growth.
        const VecXd* source = &m_rd.getGrowthFactor();
        if ((int)m_parameters["Flow type"] == 1)
            source = &m_rd.getActivator();
        if (m_stokes.solve( m_lset, *source, m_flow ))
            return -1;
        std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;


        std::cout << "* Update the interface, signaling center positions." << std::endl; tic( timeStart );
        if (m_lset.update( m_flow ))
            return -1;
        if (m_knots.updateSignalingCenters( m_lset, m_flow ))
            return -1;
        std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;


        std::cout << "* Solving morphogens." << std::endl; tic( timeStart );
        if (m_rd.solve( m_lset, m_knots ))
            return -1;
        if (m_knots.addKnots( m_lset, m_rd.getActivator() ))
            return -1;
        std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;


        std::cout << "* Writing results." << std::endl; tic( timeStart );
        if (storeModelState_(i) < 0) {
            std::cerr << "Error: Cannot write model state to output files." << std::endl;
            return -1;
        }
        std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;

        std::cout << "Iteration took " << toc( timeIter ) << " ms." << std::endl;
    }

    return 0;
}


//
//  PRIVATE METHODS
//


int Growth::storeModelState_( const int i )
{
    std::string fname = std::to_string(i) + "_" + std::to_string(m_id);
    int exportData = cfish::INTERFACE | cfish::KNOTS | cfish::BOUNDARY | cfish::INTERIOR;
    if (cfish::Write_Wavefront_obj( m_lset, m_rd, m_knots, fname, exportData ))
        return -1;

    if (cfish::Write_data( m_lset, m_rd, m_flow, fname, exportData ))
        return -1;

    return 0;
}
