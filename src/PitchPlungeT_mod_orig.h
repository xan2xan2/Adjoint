#ifndef _PITCHPLUNGET_MOD_H_
#define _PITCHPLUNGET_MOD_H_

#include <array>
#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file PitchPlunge.h
    \class PitchPlunge

    \brief Subclass of Motion, for a pitching and plunging body
    \replace frequency f with period T

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 29 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class PitchPlungeT_mod : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
    PitchPlungeT_mod(
        double dt,
        double *theta_in,
        double *y_in,
        int size
        ) :
        _dt(dt),
        _theta_in(theta_in),
        _y_in(y_in),
        _size(size) {

    }

    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {

        //////////////
        double ind_temp = std::round(time / _dt);
        int index = (int) ind_temp;
        if (index < 0){
          index = index + _size;
        }
        double y = _y_in[index];
        double ydot;
        if (index < 1){
          ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);
        }
        else if(index >= _size){
          ydot = (_y_in[index] - _y_in[index-1]) / (_dt);
        }
        else{
          ydot = (_y_in[index] - _y_in[index-1]) / (_dt);
        }

        // cout << " index ... " << index << endl;
        // cout << " size ... " << _size << endl;
        // cout << " y_minus ... " << _y_in[index-1] << endl;
        // cout << " y ... " << _y_in[index] << endl;
        // cout << " y_plus ... " << _y_in[index+1] << endl;
        // cout << " time ... " << time << endl;
        // cout << " ydot ... " << ydot << endl;


        //////////////

        double theta = _theta_in[index];
        double thetadot;
        if (index < 1){
          thetadot = (_theta_in[index + 1] - _theta_in[index]) / (_dt);
        }
        else if (index >= _size){
          thetadot = (_theta_in[index] - _theta_in[index-1]) / (_dt);
        }
        else{
          thetadot = (_theta_in[index] - _theta_in[index-1]) / (_dt);
        }

        // cout << " theta_minus ... " << _theta_in[index-1] << endl;
        // cout << " theta ... " << _theta_in[index] << endl;
        // cout << " theta_plus ... " << _theta_in[index+1] << endl;
        // cout << " thetadot ... " << thetadot << endl;
		//cout << "ydot_true..." << ydot << endl;
        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
		//                 x           x_dot   x_dot = 1 ?

    }
	///////////////////////////////////
	inline double getydot(double time) const {
    //////////////
    double ind_temp = std::round(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double ydot;
    if (index < 1){
      ydot = (_y_in[index + 1] - _y_in[index]) / (_dt);
      cout << " case1 ... " << endl;
    }
    else if (index >= _size){
      ydot = (_y_in[index] - _y_in[index-1]) / (_dt);
      cout << " case2 ... " << endl;
    }
    else{
      ydot = (_y_in[index] - _y_in[index-1]) / (_dt);
      cout << " case3 ... " << endl;
    }

    cout << " index ... " << index << endl;
    cout << " size ... " << _size << endl;
    cout << " y_minus ... " << _y_in[index-1] << endl;
    cout << " y ... " << _y_in[index] << endl;
    cout << " y_plus ... " << _y_in[index+1] << endl;
    cout << " time ... " << time << endl;
    cout << " ydot ... " << ydot << endl;
		return ydot;
	}

  inline double gety(double time) const {
    //////////////
    double ind_temp = std::round(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double y = _y_in[index];
    return y;
  }

	inline double gettheta(double time) const {
		// double theta = _pitchAmp * sin(time / _pitchT + _pitchPhase);
    //////////////
    double ind_temp = std::round(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double theta = _theta_in[index];
		return theta;
	}

	inline double getthetadot(double time) const {
		// double thetadot = _pitchAmp / _pitchT * cos(time / _pitchT + _pitchPhase);
    //////////////
    double ind_temp = std::round(time / _dt);
    int index = (int) ind_temp;
    if (index < 0){
      index = index + _size;
    }
    double thetadot;
    if (index < 1){
      thetadot = (_theta_in[index + 1] - _theta_in[index]) / (_dt);
    }
    else if (index >= _size){
      thetadot = (_theta_in[index] - _theta_in[index-1]) / (_dt);
    }
    else{
      thetadot = (_theta_in[index] - _theta_in[index-1]) / (_dt);
    }
		return thetadot;
	}

	/////////////////////////////////////
    inline Motion* clone() const {
        double twopi = 8. * atan(1.);
        return new PitchPlungeT_mod(
            _dt,
            _theta_in,
            _y_in,
            _size
        );
    };

private:
    double _dt;
    double *_theta_in;
    double *_y_in;
    int _size;

};

} // namespace ibpm

#endif /* _PITCHPLUNGET_H_ */
