#include <iostream>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#include "leveque.hpp"
#include "WriteData.hpp"
#include "fem/TetraMesh.hpp"
#include "common/Utils.hpp"

using namespace FEM;


int LeVeque::run( const int iterations )
{
    for (int i=1; i<iterations+1; i++) {
        chrono_time timeIter;
        std::cout << std::endl << "Iteration " << i << std::endl; tic( timeIter );


        chrono_time timeStart;
        std::cout << "* Computing flow field." << std::endl; tic( timeStart );
        m_flow.solve( m_lset, m_parameters, i-1, m_flowField );
        std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;


        std::cout << "* Advecting the interface." << std::endl; tic( timeStart );
        if (m_lset.update( m_flowField )) {
            std::cerr << "Unknown error updating level set." << std::endl;
            return -1;
        }
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


int LeVeque::init( const Parameters& parameters, const int id )
{
    m_parameters = parameters;
    m_id = id;

    chrono_time timeInit; 
    std::cout << std::endl << "** Initializing..." << std::endl; tic( timeInit );

    // Assuming hyperthreading is enabled. For best performance we want to use
    // the number of physical cores rather than the number of threads.
    int nThreads = std::thread::hardware_concurrency() / 2;
    std::cout << "Number of threads: " << nThreads << std::endl;
    Eigen::setNbThreads( nThreads );


    chrono_time timeStart;
    std::cout << "1/2 : Initializing level set..." << std::endl; tic( timeStart );
    if (m_lset.init( m_parameters ))
        return -1;
    std::cout << "  Mesh size: " << m_lset.getMesh().sizeInMemory() << " MB." << std::endl;
    std::cout << "  Time elapsed: " << toc( timeStart ) << " ms." << std::endl;


    std::cout << "2/3 : Initializing variables..." << std::endl; tic( timeStart );
    size_t n = m_lset.getMesh().p().size();
    m_flowField.setZero(4*n);        // 3D velocity + pressure
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



int LeVeque::storeModelState_( const int i )
{
    std::string fname = std::to_string(i) + "_" + std::to_string(m_id);

    int exportData = cfish::INTERFACE;
    if (cfish::Write_Wavefront_obj( m_lset, fname, exportData ))
        return -1;

    return 0;
}
