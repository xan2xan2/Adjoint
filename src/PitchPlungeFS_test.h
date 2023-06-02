#ifndef _PITCHPLUNGEFS_H_TEST_
#define _PITCHPLUNGEFS_H_TEST_

#include "Motion.h"
#include "TangentSE2.h"
#include <math.h>

namespace ibpm {

/*!
    \file PitchPlunge.h
    \class PitchPlunge

    \brief Subclass of Motion, for a pitching and plunging body

    \author Clancy Rowley
    \author $LastChangedBy$
    \date 29 Aug 2008
    \date $LastChangedDate$
    \version $Revision$
*/

class PitchPlungeFS_test : public Motion {
public:

    /// \brief Define a Motion corresponding to sinusoidal pitching and
    /// plunging, centered about the origin:
    ///    y(t)      = A_y      sin( 2*pi * f_y t + phi_y)
    ///    \theta(t) = A_\theta sin( 2*pi * f_theta t + phi_theta)
	PitchPlungeFS_test(
		//
		double pitchAmplitude_sin0,
		double pitchWaveNum_sin0,

		double pitchAmplitude_sin1,
		double pitchWaveNum_sin1,

		double pitchAmplitude_sin2,
		double pitchWaveNum_sin2,

		double pitchAmplitude_sin3,
		double pitchWaveNum_sin3,

		double pitchAmplitude_sin4,
		double pitchWaveNum_sin4,
		//
		double pitchAmplitude_cos0,
		double pitchWaveNum_cos0,

		double pitchAmplitude_cos1,
		double pitchWaveNum_cos1,

		double pitchAmplitude_cos2,
		double pitchWaveNum_cos2,

		double pitchAmplitude_cos3,
		double pitchWaveNum_cos3,

		double pitchAmplitude_cos4,
		double pitchWaveNum_cos4,
		//

        double pitchPhase,

		//
		double plungeAmplitude_sin0,
		double plungeWaveNum_sin0,

		double plungeAmplitude_sin1,
		double plungeWaveNum_sin1,

		double plungeAmplitude_sin2,
		double plungeWaveNum_sin2,

		double plungeAmplitude_sin3,
		double plungeWaveNum_sin3,

		double plungeAmplitude_sin4,
		double plungeWaveNum_sin4,
		//
		double plungeAmplitude_cos0,
		double plungeWaveNum_cos0,

		double plungeAmplitude_cos1,
		double plungeWaveNum_cos1,

		double plungeAmplitude_cos2,
		double plungeWaveNum_cos2,

		double plungeAmplitude_cos3,
		double plungeWaveNum_cos3,

		double plungeAmplitude_cos4,
		double plungeWaveNum_cos4,
		
		//
		double plungePhase

	) :
		_pitchAmp_sin0(pitchAmplitude_sin0),
		_pitchWN_sin0(pitchWaveNum_sin0),

		_pitchAmp_sin1(pitchAmplitude_sin1),
		_pitchWN_sin1(pitchWaveNum_sin1),

		_pitchAmp_sin2(pitchAmplitude_sin2),
		_pitchWN_sin2(pitchWaveNum_sin2),

		_pitchAmp_sin3(pitchAmplitude_sin3),
		_pitchWN_sin3(pitchWaveNum_sin3),

		_pitchAmp_sin4(pitchAmplitude_sin4),
		_pitchWN_sin4(pitchWaveNum_sin4),

		//
		_pitchAmp_cos0(pitchAmplitude_cos0),
		_pitchWN_cos0(pitchWaveNum_cos0),

		_pitchAmp_cos1(pitchAmplitude_cos1),
		_pitchWN_cos1(pitchWaveNum_cos1),

		_pitchAmp_cos2(pitchAmplitude_cos2),
		_pitchWN_cos2(pitchWaveNum_cos2),

		_pitchAmp_cos3(pitchAmplitude_cos3),
		_pitchWN_cos3(pitchWaveNum_cos3),

		_pitchAmp_cos4(pitchAmplitude_cos4),
		_pitchWN_cos4(pitchWaveNum_cos4),

        //
		_pitchPhase(pitchPhase),

        //
		_plungeAmp_sin0(plungeAmplitude_sin0),
		_plungeWN_sin0(plungeWaveNum_sin0),

		_plungeAmp_sin1(plungeAmplitude_sin1),
		_plungeWN_sin1(plungeWaveNum_sin1),

		_plungeAmp_sin2(plungeAmplitude_sin2),
		_plungeWN_sin2(plungeWaveNum_sin2),

		_plungeAmp_sin3(plungeAmplitude_sin3),
		_plungeWN_sin3(plungeWaveNum_sin3),

		_plungeAmp_sin4(plungeAmplitude_sin4),
		_plungeWN_sin4(plungeWaveNum_sin4),

       //
		_plungeAmp_cos0(plungeAmplitude_cos0),
		_plungeWN_cos0(plungeWaveNum_cos0),

		_plungeAmp_cos1(plungeAmplitude_cos1),
		_plungeWN_cos1(plungeWaveNum_cos1),

		_plungeAmp_cos2(plungeAmplitude_cos2),
		_plungeWN_cos2(plungeWaveNum_cos2),

		_plungeAmp_cos3(plungeAmplitude_cos3),
		_plungeWN_cos3(plungeWaveNum_cos3),

		_plungeAmp_cos4(plungeAmplitude_cos4),
		_plungeWN_cos4(plungeWaveNum_cos4),


        //
		_plungePhase(plungePhase)
		 {
		double twopi = 8. * atan(1.);

			_pitchWN_sin0 /= (twopi);
			_pitchWN_sin1 /= (twopi);
			_pitchWN_sin2 /= (twopi);
			_pitchWN_sin3 /= (twopi);
			_pitchWN_sin4 /= (twopi);
			//
			_pitchWN_cos0 /= (twopi);
			_pitchWN_cos1 /= (twopi);
			_pitchWN_cos2 /= (twopi);
			_pitchWN_cos3 /= (twopi);
			_pitchWN_cos4 /= (twopi);
			//
			_plungeWN_sin0 /= (twopi);
			_plungeWN_sin1 /= (twopi);
			_plungeWN_sin2 /= (twopi);
			_plungeWN_sin3 /= (twopi);
			_plungeWN_sin4 /= (twopi);
			//
			_plungeWN_cos0 /= (twopi);
			_plungeWN_cos1 /= (twopi);
			_plungeWN_cos2 /= (twopi);
			_plungeWN_cos3 /= (twopi);
			_plungeWN_cos4 /= (twopi);
	

	}
		//int size = 5;
		//
		//double A = _pitchAmp_sin0;
		//double pitchAmp_sin = [A, _pitchAmp_sin1, _pitchAmp_sin2, _pitchAmp_sin3, _pitchAmp_sin4];
		//double pitchWN_sin = [_pitchWN_sin0, _pitchWN_sin1, _pitchWN_sin2, _pitchWN_sin3, _pitchWN_sin4];
  //      double plungeAmp_sin = [_plungeAmp_sin0, _plungeAmp_sin1, _plungeAmp_sin2, _plungeAmp_sin3, _plungeAmp_sin4];
		//double plungeWN_sin = [_plungeWN_sin0, _plungeWN_sin1, _plungeWN_sin2, _plungeWN_sin3, _plungeWN_sin4];

		//double pitchAmp_cos = [_pitchAmp_cos0, _pitchAmp_cos1, _pitchAmp_cos2, _pitchAmp_cos3, _pitchAmp_cos4];
		//double pitchWN_cos = [_pitchWN_cos0, _pitchWN_cos1, _pitchWN_cos2, _pitchWN_cos3, _pitchWN_cos4];
		//double plungeAmp_cos = [_plungeAmp_cos0, _plungeAmp_cos1, _plungeAmp_cos2, _plungeAmp_cos3, _plungeAmp_cos4];
		//double plungeWN_cos = [_plungeWN_cos0, _plungeWN_cos1, _plungeWN_cos2, _plungeWN_cos3, _plungeWN_cos4];

    /// Returns transformation for sinusoidal pitch/plunge:
    ///    (0, y(t), theta(t), 0, ydot(t), thetadot(t))
    inline TangentSE2 getTransformation(double time) const {

		double theta = 0.0;
		double thetadot = 0.0;
		double y = 0.0;
		double ydot = 0.0;
		  // sin
			
				theta = _pitchAmp_sin0 * sin(1. / _pitchWN_sin0 * time + _pitchPhase) + theta;
				theta = _pitchAmp_sin1 * sin(1. / _pitchWN_sin1 * time + _pitchPhase) + theta;
				theta = _pitchAmp_sin2 * sin(1. / _pitchWN_sin2 * time + _pitchPhase) + theta;
				theta = _pitchAmp_sin3 * sin(1. / _pitchWN_sin3 * time + _pitchPhase) + theta;
				theta = _pitchAmp_sin4 * sin(1. / _pitchWN_sin4 * time + _pitchPhase) + theta;

				thetadot = _pitchAmp_sin0 / _pitchWN_sin0 * cos(1. / _pitchWN_sin0 * time + _pitchPhase) + thetadot;
				thetadot = _pitchAmp_sin1 / _pitchWN_sin1 * cos(1. / _pitchWN_sin1 * time + _pitchPhase) + thetadot;
				thetadot = _pitchAmp_sin2 / _pitchWN_sin2 * cos(1. / _pitchWN_sin2 * time + _pitchPhase) + thetadot;
				thetadot = _pitchAmp_sin3 / _pitchWN_sin3 * cos(1. / _pitchWN_sin3 * time + _pitchPhase) + thetadot;
				thetadot = _pitchAmp_sin4 / _pitchWN_sin4 * cos(1. / _pitchWN_sin4 * time + _pitchPhase) + thetadot;

				y = _plungeAmp_sin0 * sin(1. / _plungeWN_sin0 * time + _plungePhase) + y;
				y = _plungeAmp_sin1 * sin(1. / _plungeWN_sin1 * time + _plungePhase) + y;
				y = _plungeAmp_sin2 * sin(1. / _plungeWN_sin2 * time + _plungePhase) + y;
				y = _plungeAmp_sin3 * sin(1. / _plungeWN_sin3 * time + _plungePhase) + y;
				y = _plungeAmp_sin4 * sin(1. / _plungeWN_sin4 * time + _plungePhase) + y;


				ydot = _plungeAmp_sin0 / _plungeWN_sin0 * cos(1. / _plungeWN_sin0 * time + _plungePhase) + ydot;
				ydot = _plungeAmp_sin1 / _plungeWN_sin1 * cos(1. / _plungeWN_sin1 * time + _plungePhase) + ydot;
				ydot = _plungeAmp_sin2 / _plungeWN_sin2 * cos(1. / _plungeWN_sin2 * time + _plungePhase) + ydot;
				ydot = _plungeAmp_sin3 / _plungeWN_sin3 * cos(1. / _plungeWN_sin3 * time + _plungePhase) + ydot;
				ydot = _plungeAmp_sin4 / _plungeWN_sin4 * cos(1. / _plungeWN_sin4 * time + _plungePhase) + ydot;

			
		// cos
	
				theta = _pitchAmp_cos0 * cos(1. / _pitchWN_cos0 * time + _pitchPhase) + theta;
				theta = _pitchAmp_cos1 * cos(1. / _pitchWN_cos1 * time + _pitchPhase) + theta;
				theta = _pitchAmp_cos2 * cos(1. / _pitchWN_cos2 * time + _pitchPhase) + theta;
				theta = _pitchAmp_cos3 * cos(1. / _pitchWN_cos3 * time + _pitchPhase) + theta;
				theta = _pitchAmp_cos4 * cos(1. / _pitchWN_cos4 * time + _pitchPhase) + theta;


				thetadot = -_pitchAmp_cos0 / _pitchWN_cos0 * sin(1. / _pitchWN_cos0 * time + _pitchPhase) + thetadot;
				thetadot = -_pitchAmp_cos1 / _pitchWN_cos1 * sin(1. / _pitchWN_cos1 * time + _pitchPhase) + thetadot;
				thetadot = -_pitchAmp_cos2 / _pitchWN_cos2 * sin(1. / _pitchWN_cos2 * time + _pitchPhase) + thetadot;
				thetadot = -_pitchAmp_cos3 / _pitchWN_cos3 * sin(1. / _pitchWN_cos3 * time + _pitchPhase) + thetadot;
				thetadot = -_pitchAmp_cos4 / _pitchWN_cos4 * sin(1. / _pitchWN_cos4 * time + _pitchPhase) + thetadot;

				y = _plungeAmp_cos0 * cos(1. / _plungeWN_cos0 * time + _plungePhase) + y;
				y = _plungeAmp_cos1 * cos(1. / _plungeWN_cos1 * time + _plungePhase) + y;
				y = _plungeAmp_cos2 * cos(1. / _plungeWN_cos2 * time + _plungePhase) + y;
				y = _plungeAmp_cos3 * cos(1. / _plungeWN_cos3 * time + _plungePhase) + y;
				y = _plungeAmp_cos4 * cos(1. / _plungeWN_cos4 * time + _plungePhase) + y;


				ydot = -_plungeAmp_cos0 / _plungeWN_cos0 * sin(1. / _plungeWN_cos0 * time + _plungePhase) + ydot;
				ydot = -_plungeAmp_cos1 / _plungeWN_cos1 * sin(1. / _plungeWN_cos1 * time + _plungePhase) + ydot;
				ydot = -_plungeAmp_cos2 / _plungeWN_cos2 * sin(1. / _plungeWN_cos2 * time + _plungePhase) + ydot;
				ydot = -_plungeAmp_cos3 / _plungeWN_cos3 * sin(1. / _plungeWN_cos3 * time + _plungePhase) + ydot;
				ydot = -_plungeAmp_cos4 / _plungeWN_cos4 * sin(1. / _plungeWN_cos4 * time + _plungePhase) + ydot;

			

        return TangentSE2( 0, y, theta, 0, ydot, thetadot );
    }



	///////////////////////////////////
	inline double getydot(double time) const {

		double ydot = 0.0;
		//sin
		ydot = _plungeAmp_sin0 / _plungeWN_sin0 * cos(1. / _plungeWN_sin0 * time + _plungePhase) + ydot;
		ydot = _plungeAmp_sin1 / _plungeWN_sin1 * cos(1. / _plungeWN_sin1 * time + _plungePhase) + ydot;
		ydot = _plungeAmp_sin2 / _plungeWN_sin2 * cos(1. / _plungeWN_sin2 * time + _plungePhase) + ydot;
		ydot = _plungeAmp_sin3 / _plungeWN_sin3 * cos(1. / _plungeWN_sin3 * time + _plungePhase) + ydot;
		ydot = _plungeAmp_sin4 / _plungeWN_sin4 * cos(1. / _plungeWN_sin4 * time + _plungePhase) + ydot;
		//cos
		ydot = -_plungeAmp_cos0 / _plungeWN_cos0 * sin(1. / _plungeWN_cos0 * time + _plungePhase) + ydot;
		ydot = -_plungeAmp_cos1 / _plungeWN_cos1 * sin(1. / _plungeWN_cos1 * time + _plungePhase) + ydot;
		ydot = -_plungeAmp_cos2 / _plungeWN_cos2 * sin(1. / _plungeWN_cos2 * time + _plungePhase) + ydot;
		ydot = -_plungeAmp_cos3 / _plungeWN_cos3 * sin(1. / _plungeWN_cos3 * time + _plungePhase) + ydot;
		ydot = -_plungeAmp_cos4 / _plungeWN_cos4 * sin(1. / _plungeWN_cos4 * time + _plungePhase) + ydot;

		return ydot;
	}

	inline double gettheta(double time) const {
		double theta = 0.0;
		// sin
		theta = _pitchAmp_sin0 * sin(1. / _pitchWN_sin0 * time + _pitchPhase) + theta;
		theta = _pitchAmp_sin1 * sin(1. / _pitchWN_sin1 * time + _pitchPhase) + theta;
		theta = _pitchAmp_sin2 * sin(1. / _pitchWN_sin2 * time + _pitchPhase) + theta;
		theta = _pitchAmp_sin3 * sin(1. / _pitchWN_sin3 * time + _pitchPhase) + theta;
		theta = _pitchAmp_sin4 * sin(1. / _pitchWN_sin4 * time + _pitchPhase) + theta;
		// cos
		theta = _pitchAmp_cos0 * cos(1. / _pitchWN_cos0 * time + _pitchPhase) + theta;
		theta = _pitchAmp_cos1 * cos(1. / _pitchWN_cos1 * time + _pitchPhase) + theta;
		theta = _pitchAmp_cos2 * cos(1. / _pitchWN_cos2 * time + _pitchPhase) + theta;
		theta = _pitchAmp_cos3 * cos(1. / _pitchWN_cos3 * time + _pitchPhase) + theta;
		theta = _pitchAmp_cos4 * cos(1. / _pitchWN_cos4 * time + _pitchPhase) + theta;
		return theta;
	}

	inline double getthetadot(double time) const {
	
		double thetadot = 0.0;

		// sin
		thetadot = _pitchAmp_sin0 / _pitchWN_sin0 * cos(1. / _pitchWN_sin0 * time + _pitchPhase) + thetadot;
		thetadot = _pitchAmp_sin1 / _pitchWN_sin1 * cos(1. / _pitchWN_sin1 * time + _pitchPhase) + thetadot;
		thetadot = _pitchAmp_sin2 / _pitchWN_sin2 * cos(1. / _pitchWN_sin2 * time + _pitchPhase) + thetadot;
		thetadot = _pitchAmp_sin3 / _pitchWN_sin3 * cos(1. / _pitchWN_sin3 * time + _pitchPhase) + thetadot;
		thetadot = _pitchAmp_sin4 / _pitchWN_sin4 * cos(1. / _pitchWN_sin4 * time + _pitchPhase) + thetadot;
		// cos
		thetadot = -_pitchAmp_cos0 / _pitchWN_cos0 * sin(1. / _pitchWN_cos0 * time + _pitchPhase) + thetadot;
		thetadot = -_pitchAmp_cos1 / _pitchWN_cos1 * sin(1. / _pitchWN_cos1 * time + _pitchPhase) + thetadot;
		thetadot = -_pitchAmp_cos2 / _pitchWN_cos2 * sin(1. / _pitchWN_cos2 * time + _pitchPhase) + thetadot;
		thetadot = -_pitchAmp_cos3 / _pitchWN_cos3 * sin(1. / _pitchWN_cos3 * time + _pitchPhase) + thetadot;
		thetadot = -_pitchAmp_cos4 / _pitchWN_cos4 * sin(1. / _pitchWN_cos4 * time + _pitchPhase) + thetadot;
		return thetadot;
	}
	/////////////////////////////////////

    inline Motion* clone() const {
		double twopi = 8. * atan(1.);
        return new PitchPlungeFS_test(
			_pitchAmp_sin0,
			_pitchWN_sin0 * (twopi),
			
			_pitchAmp_sin1,
			_pitchWN_sin1 * (twopi),

			_pitchAmp_sin2,
			_pitchWN_sin2 * (twopi),

			_pitchAmp_sin3,
			_pitchWN_sin3 * (twopi),

			_pitchAmp_sin4,
			_pitchWN_sin4 * (twopi),
			//
			_pitchAmp_cos0,
			_pitchWN_cos0* (twopi),

			_pitchAmp_cos1,
			_pitchWN_cos1* (twopi),

			_pitchAmp_cos2,
			_pitchWN_cos2* (twopi),

			_pitchAmp_cos3,
			_pitchWN_cos3* (twopi),

			_pitchAmp_cos4,
			_pitchWN_cos4* (twopi),
	        //
			_pitchPhase,
			//
			_plungeAmp_sin0,
			_plungeWN_sin0 * (twopi),

			_plungeAmp_sin1,
			_plungeWN_sin1* (twopi),

			_plungeAmp_sin2,
			_plungeWN_sin2* (twopi),

			_plungeAmp_sin3,
			_plungeWN_sin3* (twopi),

			_plungeAmp_sin4,
			_plungeWN_sin4* (twopi),
			//
			_plungeAmp_cos0,
			_plungeWN_cos0* (twopi),

			_plungeAmp_cos1,
			_plungeWN_cos1* (twopi),

			_plungeAmp_cos2,
			_plungeWN_cos2* (twopi),

			_plungeAmp_cos3,
			_plungeWN_cos3* (twopi),

			_plungeAmp_cos4,
			_plungeWN_cos4* (twopi),
			//

			_plungePhase

        );
    };

private:
    double _pitchAmp_sin0;
    double _pitchWN_sin0;

	double _pitchAmp_sin1;
	double _pitchWN_sin1;

	double _pitchAmp_sin2;
	double _pitchWN_sin2;

	double _pitchAmp_sin3;
	double _pitchWN_sin3;

	double _pitchAmp_sin4;
	double _pitchWN_sin4;
	//

	double _pitchAmp_cos0;
	double _pitchWN_cos0;

	double _pitchAmp_cos1;
	double _pitchWN_cos1;

	double _pitchAmp_cos2;
	double _pitchWN_cos2;

	double _pitchAmp_cos3;
	double _pitchWN_cos3;

	double _pitchAmp_cos4;
	double _pitchWN_cos4;

	//
    double _pitchPhase;
	//
    double _plungeAmp_sin0;
    double _plungeWN_sin0;

	double _plungeAmp_sin1;
	double _plungeWN_sin1;

	double _plungeAmp_sin2;
	double _plungeWN_sin2;

	double _plungeAmp_sin3;
	double _plungeWN_sin3;

	double _plungeAmp_sin4;
	double _plungeWN_sin4;
	//
	double _plungeAmp_cos0;
	double _plungeWN_cos0;

	double _plungeAmp_cos1;
	double _plungeWN_cos1;

	double _plungeAmp_cos2;
	double _plungeWN_cos2;

	double _plungeAmp_cos3;
	double _plungeWN_cos3;

	double _plungeAmp_cos4;
	double _plungeWN_cos4;

	//
    double _plungePhase;

};

} // namespace ibpm

#endif /* _PITCHPLUNGE_H_ */
