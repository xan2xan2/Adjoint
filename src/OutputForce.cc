// OutputForce.cc
//
// Description:
// Implementation of output routine for writing a list of force coefficients.
//
// Author(s):
// Clancy Rowley
// Steve Brunton
//
// Date: 22 Aug 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "OutputForce.h"
#include "BaseFlow.h"
#include "State.h"
#include "Output.h"
#include "VectorOperations.h"
#include <stdio.h>
#include <string>

//using namespace std;
using std::string;

namespace ibpm {

OutputForce::OutputForce(string filename) :
    _filename( filename )
{}

bool OutputForce::init() {
    _fp = fopen( _filename.c_str(), "w" );
    if ( _fp == NULL ) return false;
    else return true;
}

bool OutputForce::cleanup() {
    bool status = true;
    if ( _fp != NULL ) {
        status = fclose( _fp );
    }
    return status;
}

// Method to compute lift, drag from state (x), angle of attack (alpha), and freestream velocity (mag)
bool OutputForce::doOutput( const double alpha, const double theta, const double mag, const State& x, const double zM) {
    double xF, yF;      // force in x and y direction in domain
    double drag, lift;  // actual drag and lift, wrt alpha
    x.computeNetForce( xF, yF );

    ////////////////////////
    // double ydot = motion.getydot(x.time);
    // double theta = motion.gettheta(x.time);
    // double thetadot = motion.getthetadot(x.time);

    ////////////////////

    drag = xF * cos(theta) - yF * sin(theta);
    lift = xF * sin(theta) + yF * cos(theta);
    ////////////////////////
    // drag = xF * cos(alpha) + yF * sin(alpha);
    // lift = xF * -1.*sin(alpha) + yF * cos(alpha);

    // Convert forces to lift and drag coefficients:
    // If L_d is dimensional lift, then in the nondimensionalization of the
    // code (lengths by c, density by rho, velocity by U), we have
    //    L = L_d / (c rho U^2)
    // so
    //    C_L = L_d / (1/2 rho U^2)
    //        = 2 L
    // drag *= 2.;
    // lift *= 2.;

    if ( _fp == NULL ) return false;
    fprintf( _fp, "%5d %.5e %.5e %.5e %.5e %.5e\n", x.timestep, x.time, theta, drag, lift, zM );
    fflush( _fp );

    return true;
}

// If no other information is provided, assume zero angle of attack, unity freestrem velocity
bool OutputForce::doOutput(const State& x) {
    double alpha = 0.;
    double mag = 1.;
    double theta = 0.;
    double zM = 0.;

    return doOutput(alpha, theta, mag, x, zM);
}

bool OutputForce::doOutput(const BaseFlow& q, const State& x) {
    double alpha = q.getAlpha();
    double mag = q.getMag();
    // double theta = q.getTheta();
    double theta = 0.;
    double zM = 0.;

    return doOutput(alpha, theta, mag, x, zM);
}
////////
bool OutputForce::doOutput(const double theta, const BaseFlow& q, const State& x, const double zM) {
    double alpha = q.getAlpha();
    double mag = q.getMag();
    // double theta = motion.gettheta(x.time);
    // double theta = q.getTheta();

    return doOutput(alpha, theta, mag, x, zM);
}
////////
} // namespace ibpm
