/* -------------------------------------------------------------------------- *
 *                   OpenSim:  StaticOptimization_GHTarget.cpp                   *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2012 Stanford University and the Authors                *
 * Author(s): Frank C. Anderson                                               *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/* Note: This code was originally developed by Realistic Dynamics Inc. 
 * Author: Frank C. Anderson 
 */


//=============================================================================
// INCLUDES
//=============================================================================
#include <stdlib.h>
#include <stdio.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/ActivationFiberLengthMuscle.h>
#include <OpenSim/Simulation/Model/ForceSet.h>
#include <OpenSim/Simulation/Model/BodySet.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>
#include "StaticOptimization_GHTarget.h"
#include <iostream>

using namespace OpenSim;
using namespace std;
using SimTK::Vector;
using SimTK::Matrix;
using SimTK::Real;

#define USE_LINEAR_CONSTRAINT_MATRIX
#define USE_LINEAR_GHFORCE_MATRIX

const double StaticOptimization_GHTarget::SMALLDX = 1.0e-14;
//const double StaticOptimization_GHTarget::_activationExponent = 2.0;
 
//==============================================================================
// CONSTRUCTOR
//==============================================================================
//______________________________________________________________________________
/**
 * Construct an optimization target.
 *
 * @param aNP The number of parameters.
 * @param aNC The number of constraints.
 * @param aT Current time in the integration.
 * @param aX Current control values.
 * @param aY Current states.
 * @param aDYDT Current state derivatives.
 */
StaticOptimization_GHTarget::
StaticOptimization_GHTarget(const SimTK::State& s, Model *aModel,int aNP,int aNC, bool useMusclePhysiology)
{

	// ALLOCATE STATE ARRAYS
	_recipAreaSquared.setSize(aNP);
	_recipOptForceSquared.setSize(aNP);
	_optimalForce.setSize(aNP);
	_useMusclePhysiology=useMusclePhysiology;

	setModel(*aModel);
	setNumParams(aNP);
	setNumConstraints(aNC);
	setActivationExponent(2.0);
	computeActuatorAreas(s);

	// Gather indices into speed set corresponding to the unconstrained degrees of freedom (for which we will set acceleration constraints)
	_accelerationIndices.setSize(0);
	const CoordinateSet& coordSet = _model->getCoordinateSet();
	for(int i=0; i<coordSet.getSize(); i++) {
		const Coordinate& coord = coordSet.get(i);
		if(!coord.isConstrained(s)) {
			_accelerationIndices.append(i);
		}
	}

}


//==============================================================================
// CONSTRUCTION
//==============================================================================
bool StaticOptimization_GHTarget::
prepareToOptimize(SimTK::State& s, double *x)
{
	// Keep around a "copy" of the state so we can use it in objective function 
	// in cases where we're tracking states
	_currentState = &s;
	
	// COMPUTE MAX ISOMETRIC FORCE
	const ForceSet& fSet = _model->getForceSet();
    
	for(int i=0, j=0;i<fSet.getSize();i++) {
 		 Actuator* act = dynamic_cast<Actuator*>(&fSet.get(i));
         if( act ) {
             double fOpt;
             Muscle *mus = dynamic_cast<Muscle*>(&fSet.get(i));
             if( mus ) {
				//ActivationFiberLengthMuscle *aflmus = dynamic_cast<ActivationFiberLengthMuscle*>(mus);
				if(mus && _useMusclePhysiology) {
					_model->setAllControllersEnabled(true);
    				fOpt = mus->calcInextensibleTendonActiveFiberForce(s, 1.0);
					_model->setAllControllersEnabled(false);
    			} else {
    				fOpt = mus->getMaxIsometricForce();
                }
             } else {
                  fOpt = act->getOptimalForce();
             }
		    _optimalForce[j++] = fOpt;
		 }
	}

#ifdef USE_LINEAR_CONSTRAINT_MATRIX
	//cout<<"Computing linear constraint matrix..."<<endl;
	int np = getNumParameters();
	int nc = getNumConstraints();

	_constraintMatrix.resize(nc,np);
	_constraintVector.resize(nc);

	Vector pVector(np), cVector(nc);

	// Build linear constraint matrix and constant constraint vector
	pVector = 0;
	computeConstraintVector(s, pVector,_constraintVector);

	for(int p=0; p<np; p++) {
		pVector[p] = 1;
		computeConstraintVector(s, pVector, cVector);
		for(int c=0; c<nc; c++) _constraintMatrix(c,p) = (cVector[c] - _constraintVector[c]);
		pVector[p] = 0;
	}
#endif

#ifdef USE_LINEAR_GHFORCE_MATRIX
	//cout<<"Computing linear constraint matrix..."<<endl;
	int nnp = getNumParameters();

	_ghforceMatrix.resize(3, nnp);
	_ghforceVector.resize(3);

	Vector pghfVector(nnp), ghfVector(3);

	// Build linear GH force matrix and constant GH fore vector
	pghfVector = 0;
	computeGHForceVector(s, pghfVector, _ghforceVector);

	for (int p = 0; p<nnp; p++) {
		pghfVector[p] = 1;
		computeGHForceVector(s, pghfVector, ghfVector);
		for (int c = 0; c<3; c++) _ghforceMatrix(c, p) = (ghfVector[c] - _ghforceVector[c]);
		pghfVector[p] = 0;
	}
#endif
	// return false to indicate that we still need to proceed with optimization
	return false;
}
//==============================================================================
// SET AND GET
//==============================================================================
//------------------------------------------------------------------------------
// MODEL
//------------------------------------------------------------------------------
///______________________________________________________________________________
/**
 * Set the model.
 *
 * @param aModel Model.
 */
void StaticOptimization_GHTarget::
setModel(Model& aModel)
{
	_model = &aModel;
}
//------------------------------------------------------------------------------
// STATES STORAGE
//------------------------------------------------------------------------------
///______________________________________________________________________________
/**
 * Set the states storage.
 *
 * @param aStatesStore States storage.
 */
void StaticOptimization_GHTarget::
setStatesStore(const Storage *aStatesStore)
{
	_statesStore = aStatesStore;
}
//------------------------------------------------------------------------------
// STATES SPLINE SET
//------------------------------------------------------------------------------
///______________________________________________________________________________
/**
 * Set the states spline set.
 *
 * @param aStatesSplineSet States spline set.
 */
void StaticOptimization_GHTarget::
setStatesSplineSet(GCVSplineSet aStatesSplineSet)
{
	_statesSplineSet = aStatesSplineSet;
}

//------------------------------------------------------------------------------
// CONTROLS
//------------------------------------------------------------------------------
///______________________________________________________________________________
/**
 * Set the number of paramters.
 *
 * The number of parameters can be set at any time.  However, the perturbation
 * sizes for the parameters (i.e., _dx) is destroyed.  Therefore, the
 * perturbation sizes must be reset.
 *
 * @param aNP Number of parameters.
 * @see setDX()
 */
void StaticOptimization_GHTarget::
setNumParams(const int aNP)
{
	setNumParameters(aNP);
	_dx.setSize(getNumParameters());
}

//------------------------------------------------------------------------------
// CONSTRAINTS
//------------------------------------------------------------------------------
///______________________________________________________________________________
/**
 * Set the number of constraints.
 *
 * @param aNC Number of constraints.
 */
void StaticOptimization_GHTarget::
setNumConstraints(const int aNC)
{
	// There are only linear equality constraints.
	setNumEqualityConstraints(aNC);
	setNumLinearEqualityConstraints(aNC);
}	

//------------------------------------------------------------------------------
// DERIVATIVE PERTURBATION SIZES
//------------------------------------------------------------------------------
//______________________________________________________________________________
/**
 * Set the derivative perturbation size.
 */
void StaticOptimization_GHTarget::
setDX(int aIndex,double aValue)
{
	// VALIDATE VALUE
	validatePerturbationSize(aValue);

	// SET VALUE (use get to do bounds checking)
	_dx.updElt(aIndex) = aValue;
}
//______________________________________________________________________________
/**
 * Set the derivative perturbation size for all controls.
 */
void StaticOptimization_GHTarget::
setDX(double aValue)
{
	// VALIDATE VALUE
	validatePerturbationSize(aValue);

	// SET VALUE
	for(int i=0;i<getNumParameters();i++) _dx.updElt(i) = aValue;
}
//______________________________________________________________________________
/**
 * Get the derivative perturbation size.
 */
double StaticOptimization_GHTarget::
getDX(int aIndex)
{
	return _dx.get(aIndex);
}
//______________________________________________________________________________
/**
 * Get a pointer to the vector of derivative perturbation sizes.
 */
double* StaticOptimization_GHTarget::
getDXArray()
{
	return &_dx[0];
}

//______________________________________________________________________________
/**
 * Get an optimal force.
 */
void StaticOptimization_GHTarget::
getActuation(SimTK::State& s, const SimTK::Vector &parameters, SimTK::Vector &forces)
{
	//return(_optimalForce[aIndex]);
	const ForceSet& fs = _model->getForceSet();
	SimTK::Vector tempAccel(getNumConstraints());
	computeAcceleration(s, parameters, tempAccel);
	for(int i=0,j=0;i<fs.getSize();i++) {
        Actuator* act = dynamic_cast<Actuator*>(&fs.get(i));
		if( act )forces(j++) = act->getForce(s);
	}
}
//==============================================================================
// UTILITY
//==============================================================================
//______________________________________________________________________________
/**
 * Ensure that a derivative perturbation is a valid size
 */
void StaticOptimization_GHTarget::
validatePerturbationSize(double &aSize)
{
	if(aSize<SMALLDX) {
		printf("StaticOptimization_GHTarget.validatePerturbationSize: WARNING- ");
		printf("dx size too small (%le).\n",aSize);
		printf("\tResetting dx=%le.\n",SMALLDX);
		aSize = SMALLDX;
	}
}
//______________________________________________________________________________
/**
 */
void StaticOptimization_GHTarget::
printPerformance(SimTK::State& s, double *parameters)
{
	double p;
	setCurrentState( &s );
	objectiveFunc(SimTK::Vector(getNumParameters(),parameters,true),true,p);
	SimTK::Vector constraints(getNumConstraints());
	constraintFunc(SimTK::Vector(getNumParameters(),parameters,true),true,constraints);
	cout << endl;
	cout << "time = " << s.getTime() <<" Performance =" << p << 
	" Constraint violation = " << sqrt(~constraints*constraints) << endl;
}

//______________________________________________________________________________
/**
 */
void StaticOptimization_GHTarget::
computeActuatorAreas(const SimTK::State& s )
{
	// COMPUTE ACTUATOR AREAS
	ForceSet& forceSet = _model->updForceSet();
	for(int i=0, j=0;i<forceSet.getSize();i++) {
        Actuator *act = dynamic_cast<Actuator*>(&forceSet.get(i));
        if( act ) {
 		     act->setForce(s, 1.0);
    		 _recipAreaSquared[j] = act->getStress(s);
    		 _recipAreaSquared[j] *= _recipAreaSquared[j];
             j++;
        }
	}
}

//=============================================================================
// STATIC DERIVATIVES
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute derivatives of a constraint with respect to the
 * controls by central differences.
 *
 * @param dx An array of control perturbation values.
 * @param x Values of the controls at time t.
 * @param ic Index of the constraint.
 * @param dcdx The derivatives of the constraints.
 *
 * @return -1 if an error is encountered, 0 otherwize.
 */
int StaticOptimization_GHTarget::
CentralDifferencesConstraint(const StaticOptimization_GHTarget *aTarget,
	double *dx,const Vector &x,Matrix &jacobian)
{
	if(aTarget==NULL) return(-1);

	// INITIALIZE CONTROLS
	int nx = aTarget->getNumParameters(); if(nx<=0) return(-1);
	int nc = aTarget->getNumConstraints(); if(nc<=0) return(-1);
	Vector xp=x;
	Vector cf(nc),cb(nc);

	// INITIALIZE STATUS
	int status = -1;

	// LOOP OVER CONTROLS
	for(int i=0;i<nx;i++) {

		// PERTURB FORWARD
		xp[i] = x[i] + dx[i];
		status = aTarget->constraintFunc(xp,true,cf);
		if(status<0) return(status);

		// PERTURB BACKWARD
		xp[i] = x[i] - dx[i];
		status = aTarget->constraintFunc(xp,true,cb);
		if(status<0) return(status);

		// DERIVATIVES OF CONSTRAINTS
		double rdx = 0.5 / dx[i];
		for(int j=0;j<nc;j++) jacobian(j,i) = rdx*(cf[j]-cb[j]);

		// RESTORE CONTROLS
		xp[i] = x[i];
	}

	return(status);
}
//_____________________________________________________________________________
/**
 * Compute derivatives of performance with respect to the
 * controls by central differences.  Note that the gradient array should
 * be allocated as dpdx[nx].
 *
 * @param dx An array of control perturbation values.
 * @param x Values of the controls at time t.
 * @param dpdx The derivatives of the performance criterion.
 *
 * @return -1 if an error is encountered, 0 otherwize.
 */
int StaticOptimization_GHTarget::
CentralDifferences(const StaticOptimization_GHTarget *aTarget,
	double *dx,const Vector &x,Vector &dpdx)
{
	if(aTarget==NULL) return(-1);

	// CONTROLS
	int nx = aTarget->getNumParameters();  if(nx<=0) return(-1);
	Vector xp=x;

	// PERFORMANCE
	double pf,pb;

	// INITIALIZE STATUS
	int status = -1;

	// LOOP OVER CONTROLS
	for(int i=0;i<nx;i++) {

		// PERTURB FORWARD
		xp[i] = x[i] + dx[i];
		status = aTarget->objectiveFunc(xp,true,pf);
		if(status<0) return(status);

		// PERTURB BACKWARD
		xp[i] = x[i] - dx[i];
		status = aTarget->objectiveFunc(xp,true,pb);
		if(status<0) return(status);

		// DERIVATIVES OF PERFORMANCE
		double rdx = 0.5 / dx[i];
		dpdx[i] = rdx*(pf-pb);

		// RESTORE CONTROLS
		xp[i] = x[i];
	}

	return(status);
}

//==============================================================================
// PERFORMANCE AND CONSTRAINTS
//==============================================================================
//------------------------------------------------------------------------------
// PERFORMANCE
//------------------------------------------------------------------------------
//______________________________________________________________________________
/**
 * Compute performance given parameters.
 *
 * @param parameters Vector of optimization parameters.
 * @param performance Value of the performance criterion.
 * @return Status (normal termination = 0, error < 0).
 */
int StaticOptimization_GHTarget::
objectiveFunc(const Vector &parameters, const bool new_parameters, Real &performance) const
{
	//LARGE_INTEGER start;
	//LARGE_INTEGER stop;
	//LARGE_INTEGER frequency;

	//QueryPerformanceFrequency(&frequency);
	//QueryPerformanceCounter(&start);

	int na = _model->getActuators().getSize();
	double p = 0.0;
	for(int i=0;i<na;i++) {
		p +=  pow(fabs(parameters[i]),_activationExponent);
	}
//	performance = p;

#ifndef USE_LINEAR_GHFORCE_MATRIX

	SimTK::Vector force;
	computeGHForceVector(*_currentState, parameters, force);

#else

	// Use precomputed GHforce matrix
	 SimTK::Vector force = _ghforceMatrix * parameters + _ghforceVector;

#endif

	/* Take glenoid orientation into account */
	double Rgt_x[3] = { 0.9511, -0.2724, -0.1455 };
	double Rgt_y[3] = { 0.2693, 0.9622, -0.0412 };
	double Rgt_z[3] = { 0.1513, 0, 0.9885 };

	double gforce[3] = { 0 };
	/* The glenoid orientation matrix came from the original Delft model 
	 * which had x pointing laterally, y superiorly and z posteriorly */
	gforce[0] = Rgt_x[0] * force[2] + Rgt_x[1] * force[1] - Rgt_x[2] * force[0];
	gforce[1] = Rgt_y[0] * force[2] + Rgt_y[1] * force[1] - Rgt_y[2] * force[0];
	gforce[2] = Rgt_z[0] * force[2] + Rgt_z[1] * force[1] - Rgt_z[2] * force[0];

	double phi, theta;
	double p2 = 0;

	// normalize the reaction force
	double fvlength = sqrt(gforce[0] * gforce[0] + gforce[1] * gforce[1] + gforce[2] * gforce[2]);
	gforce[0] = gforce[0] / fvlength;
	gforce[1] = gforce[1] / fvlength;
	gforce[2] = gforce[2] / fvlength;

	// calculate position of glenoid cavity
	theta = asin(-gforce[1]);
	if (sqrt(gforce[0] * gforce[0] + gforce[2] * gforce[2]) == 0)
		phi = 0.0;
	else
		phi = asin(gforce[2] / sqrt(gforce[0] * gforce[0] + gforce[2] * gforce[2]));

	p2 = pow((theta / _atheta), 2) + pow((phi / _aphi), 2); // this needs to be less than 1

	double fp2 = 100*p2*((p2 - 1) + sqrt(pow((p2 - 1), 2) + _GHstabTransitionParam));
//	double fp2 = 10 * p2*p2;

	performance = p + fp2;

	return(0);
}
//______________________________________________________________________________
/**
* Compute GH force vector given parameters.
*/
void StaticOptimization_GHTarget::
computeGHForceVector(SimTK::State& s, const Vector &parameters, Vector &ghforce) const
{

	// Compute actual accelerations
	Vector actualAcceleration(getNumConstraints());
	computeAcceleration(s, parameters, actualAcceleration);

	int numBodies = _model->getNumBodies();
	/** define 2 variable length vectors of Vec3 vectors to contain calculated
	*   forces and moments for all the bodies in the model */
	SimTK::Vector_<SimTK::Vec3> allForcesVec(numBodies);
	SimTK::Vector_<SimTK::Vec3> allMomentsVec(numBodies);

	//// BodySet and JointSet and ground body index
	const BodySet& bodySet = _model->getBodySet();
	const JointSet& jointSet = _model->getJointSet();
	Body &ground = _model->getSimbodyEngine().getGroundBody();

	/* Calculate All joint reaction forces anSimTK::d moments.
	*  Applied to child bodies, expressed in ground frame.
	*  computeReactions realizes to the acceleration stage internally
	*  so you don't have to call realize in this analysis.*/
	_model->getSimbodyEngine().computeReactions(*_currentState, allForcesVec, allMomentsVec);

	/* retrieved desired joint reactions, convert to desired bodies, and convert
	*  to desired reference frames*/
	const Joint& joint = jointSet.get("glenohumeral");
	SimTK::Vec3 force = allForcesVec[jointSet.getIndex("glenohumeral")];
	// convert to parent frame
	force = -force;
	const Body& expressedInBody = bodySet.get("scapula");
	/* Express force in the scapula reference frame*/
	_model->getSimbodyEngine().transform(*_currentState, ground, force, expressedInBody, force);
	for (int c = 0; c < 3; c++) ghforce[c] = force[c];

}

//______________________________________________________________________________
/**
 * Compute the gradient of performance given parameters.
 *
 * @param parameters Vector of optimization parameters.
 * @param gradient Derivatives of performance with respect to the parameters.
 * @return Status (normal termination = 0, error < 0).
 */
int StaticOptimization_GHTarget::
gradientFunc(const Vector &parameters, const bool new_parameters, Vector &gradient) const
{
	//LARGE_INTEGER start;
	//LARGE_INTEGER stop;
	//LARGE_INTEGER frequency;

	//QueryPerformanceFrequency(&frequency);
	//QueryPerformanceCounter(&start);

//	int na = _model->getActuators().getSize();
//	for(int i=0;i<na;i++) {
//		if(parameters[i] < 0) {
//			gradient[i] =  -1.0 * _activationExponent * pow(fabs(parameters[i]),_activationExponent-1.0);
//		} else {
//			gradient[i] =  _activationExponent * pow(fabs(parameters[i]),_activationExponent-1.0);
//	}
//	}

	//QueryPerformanceCounter(&stop);
	//double duration = (double)(stop.QuadPart-start.QuadPart)/(double)frequency.QuadPart;
	//std::cout << "gradientFunc time = " << (duration*1.0e3) << " milliseconds" << std::endl;

	// 0.02 ms

	// Compute gradient 
	StaticOptimization_GHTarget::CentralDifferences(this, &_dx[0], parameters, gradient);

	return(0);
}

//------------------------------------------------------------------------------
// CONSTRAINT
//------------------------------------------------------------------------------
//______________________________________________________________________________
/**
 * Compute acceleration constraints given parameters.
 *
 * @param parameters Vector of optimization parameters.
 * @param constraints Vector of optimization constraints.
 * @return Status (normal termination = 0, error < 0).
 */
int StaticOptimization_GHTarget::
constraintFunc(const SimTK::Vector &parameters, const bool new_parameters, SimTK::Vector &constraints) const
{
	//LARGE_INTEGER start;
	//LARGE_INTEGER stop;
	//LARGE_INTEGER frequency;

	//QueryPerformanceFrequency(&frequency);
	//QueryPerformanceCounter(&start);

#ifndef USE_LINEAR_CONSTRAINT_MATRIX

	// Evaluate constraint function for all constraints and pick the appropriate component
	computeConstraintVector(parameters,constraints);

#else

	// Use precomputed constraint matrix
	//cout<<"Computing constraints assuming linear dependence..."<<endl;
	constraints = _constraintMatrix * parameters + _constraintVector;

#endif

	//QueryPerformanceCounter(&stop);
	//double duration = (double)(stop.QuadPart-start.QuadPart)/(double)frequency.QuadPart;
	//std::cout << "constraintFunc time = " << (duration*1.0e3) << " milliseconds" << std::endl;

	// 0.11 ms

	return(0);
}

//______________________________________________________________________________
/**
 * Compute all constraints given parameters.
 */
void StaticOptimization_GHTarget::
computeConstraintVector(SimTK::State& s, const Vector &parameters,Vector &constraints) const
{
	//LARGE_INTEGER start;
	//LARGE_INTEGER stop;
	//LARGE_INTEGER frequency;

	//QueryPerformanceFrequency(&frequency);
	//QueryPerformanceCounter(&start);

	// Compute actual accelerations
	Vector actualAcceleration(getNumConstraints());
	computeAcceleration(s, parameters, actualAcceleration);

	// CONSTRAINTS
	for(int i=0; i<getNumConstraints(); i++) {
		Coordinate& coord = _model->getCoordinateSet().get(_accelerationIndices[i]);
		Function& presribedFunc = _statesSplineSet.get(_statesStore->getStateIndex(coord.getSpeedName(),0));
		std::vector<int> derivComponents(1,0); //take first derivative
		double targetAcceleration = presribedFunc.calcDerivative(derivComponents,SimTK::Vector(1,s.getTime()));
		//std::cout << "computeConstraintVector:" << targetAcceleration << " - " <<  actualAcceleration[i] << endl;
		constraints[i] = targetAcceleration - actualAcceleration[i];
	}

	//QueryPerformanceCounter(&stop);
	//double duration = (double)(stop.QuadPart-start.QuadPart)/(double)frequency.QuadPart;
	//std::cout << "computeConstraintVector time = " << (duration*1.0e3) << " milliseconds" << std::endl;

	// 1.5 ms
}
//______________________________________________________________________________
/**
 * Compute the gradient of constraint given parameters.
 *
 * @param parameters Vector of parameters.
 * @param jac Derivative of constraint with respect to the parameters.
 * @return Status (normal termination = 0, error < 0).
 */
int StaticOptimization_GHTarget::
constraintJacobian(const SimTK::Vector &parameters, const bool new_parameters, SimTK::Matrix &jac) const
{
	//LARGE_INTEGER start;
	//LARGE_INTEGER stop;
	//LARGE_INTEGER frequency;

	//QueryPerformanceFrequency(&frequency);
	//QueryPerformanceCounter(&start);

#ifndef USE_LINEAR_CONSTRAINT_MATRIX

	// Compute gradient 
	StaticOptimization_GHTarget::CentralDifferencesConstraint(this,&_dx[0],parameters,jac);

#else

	// Use precomputed constraint matrix (works if constraint is linear)
	//cout<<"Computing constraint gradient assuming linear dependence..."<<endl;
	jac = _constraintMatrix;

#endif

	//QueryPerformanceCounter(&stop);
	//double duration = (double)(stop.QuadPart-start.QuadPart)/(double)frequency.QuadPart;
	//std::cout << "constraintJacobian time = " << (duration*1.0e3) << " milliseconds" << std::endl;

	// 0.01 ms

	return 0;
}
//=============================================================================
// ACCELERATION
//=============================================================================
//
void StaticOptimization_GHTarget::
computeAcceleration(SimTK::State& s, const SimTK::Vector &parameters,SimTK::Vector &rAccel) const
{
	//LARGE_INTEGER start;
	//LARGE_INTEGER stop;
	//LARGE_INTEGER frequency;

	//QueryPerformanceFrequency(&frequency);
	//QueryPerformanceCounter(&start);

	// SimTK requires that time be >= 0 when setting Discreate variables (overrideForce)
	// JACKM: Need to talk to sherm if this restriction can be removed
	double time = s.getTime();
	

    const Set<Actuator>& as = _model->getForceSet().getActuators();
	for(int i=0;i<as.getSize();i++)  {
        as.get(i).setOverrideForce(s, parameters[i] * _optimalForce[i]);
    }

	_model->getMultibodySystem().realize(s,SimTK::Stage::Acceleration);

	SimTK::Vector udot = _model->getMatterSubsystem().getUDot(s);

	for(int i=0; i<_accelerationIndices.getSize(); i++) 
		rAccel[i] = udot[_accelerationIndices[i]];

	//QueryPerformanceCounter(&stop);
	//double duration = (double)(stop.QuadPart-start.QuadPart)/(double)frequency.QuadPart;
	//std::cout << "computeAcceleration time = " << (duration*1.0e3) << " milliseconds" << std::endl;

	// 1.45 ms
}
