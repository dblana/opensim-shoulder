/* -------------------------------------------------------------------------- *
 *                        OpenSim:  JointReaction_GH.cpp                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2012 Stanford University and the Authors                *
 * Author(s): Matt S. DeMers, Ajay Seth                                       *
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


//=============================================================================
// INCLUDES
//=============================================================================
#include <iostream>
#include <string>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Actuator.h>
#include <OpenSim/Simulation/Model/BodySet.h>
#include <OpenSim/Simulation/SimbodyEngine/SimbodyEngine.h>
#include "JointReaction_GH.h"

using namespace OpenSim;
using namespace std;
using namespace SimTK;


//=============================================================================
// CONSTANTS
//=============================================================================



//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Destructor.
 * Delete any variables allocated using the "new" operator.  You will not
 * necessarily have any of these.
 */
JointReaction_GH::~JointReaction_GH()
{

}
//_____________________________________________________________________________
/**
 * Construct a JointReaction_GH instance.
 *
 * @param aModel Model for which the analysis is to be run.
 */
JointReaction_GH::JointReaction_GH(Model *aModel) :
	Analysis(aModel),
	_forcesFileName(_forcesFileNameProp.getValueStr()),
	_jointNames(_jointNamesProp.getValueStrArray()),
	_onBody(_onBodyProp.getValueStrArray()),
	_inFrame(_inFrameProp.getValueStrArray()),
	_atheta(_athetaProp.getValueDbl()),
	_aphi(_aphiProp.getValueDbl())
{
	setNull();
}
//_____________________________________________________________________________
/**
 * Construct an object from file.
 *
 * The object is constructed from the root element of the XML document.
 * The type of object is the tag name of the XML root element.
 *
 * @param aFileName File name of the document.
 */
JointReaction_GH::JointReaction_GH(const std::string &aFileName):
	Analysis(aFileName, false),
	_forcesFileName(_forcesFileNameProp.getValueStr()),
	_jointNames(_jointNamesProp.getValueStrArray()),
	_onBody(_onBodyProp.getValueStrArray()),
	_inFrame(_inFrameProp.getValueStrArray()),
	_atheta(_athetaProp.getValueDbl()),
	_aphi(_aphiProp.getValueDbl())
{
	setNull();

	// Serialize from XML
	updateFromXMLDocument();

	/* The rest will be done by setModel().
	// CONSTRUCT DESCRIPTION AND LABELS
	constructDescription();
	updateBodiesToRecord();
	constructColumnLabels();

	// STORAGE
	allocateStorage();
	*/
}

// Copy constrctor and virtual copy 
//_____________________________________________________________________________
/**
 * Copy constructor.
 *
 */
JointReaction_GH::JointReaction_GH(const JointReaction_GH &aJointReaction_GH):
	Analysis(aJointReaction_GH),
	_forcesFileName(_forcesFileNameProp.getValueStr()),
	_jointNames(_jointNamesProp.getValueStrArray()),
	_onBody(_onBodyProp.getValueStrArray()),
	_inFrame(_inFrameProp.getValueStrArray()),
	_atheta(_athetaProp.getValueDbl()),
	_aphi(_aphiProp.getValueDbl())
{
	setNull();
	// COPY TYPE AND NAME
	*this = aJointReaction_GH;
}

//=============================================================================
// OPERATORS
//=============================================================================
//-----------------------------------------------------------------------------
// ASSIGNMENT
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
/**
 * Assign this object to the values of another.
 *
 * @return Reference to this object.
 */
JointReaction_GH& JointReaction_GH::
operator=(const JointReaction_GH &aJointReaction_GH)
{
	// Base Class
	Analysis::operator=(aJointReaction_GH);

	// Member Variables
	_forcesFileName = aJointReaction_GH._forcesFileName;
	_jointNames = aJointReaction_GH._jointNames;
	_onBody.setSize(1);
	_onBody[0] = "child";
	_inFrame.setSize(1);
	_inFrame[0] = "parent";
	_atheta = aJointReaction_GH._atheta;
	_aphi = aJointReaction_GH._aphi;
	_useForceStorage = aJointReaction_GH._useForceStorage;
	_storeActuation = NULL;
	return(*this);
}

//_____________________________________________________________________________
/**
 * SetNull().
 */
void JointReaction_GH::
setNull()
{
	setAuthors("Dimitra Blana");
	setupProperties();

	// Property Default Values that are set if the associated fields are
	// omitted from the setup file.
	_forcesFileName = "";
	_useForceStorage = false;
	_jointNames.setSize(1);
	_jointNames[0] = "glenohumeral";
	_atheta = 0.7744;
	_aphi = 0.6728;
	_storeActuation = NULL;

}
//_____________________________________________________________________________
/**
 * Set up the properties for the analysis.
 *
 * The name give to each property is the tag that will be used in the XML
 * file.  The comment will appear before the property in the XML file.
 * In addition, the comments are used for tool tips in the OpenSim GUI.
 *
 * All properties are added to the property set.  Once added, they can be
 * read in and written to file.
 */
void JointReaction_GH::
setupProperties()
{

	_forcesFileNameProp.setName("forces_file");
	_forcesFileNameProp.setComment("The name of a file containing forces storage."
		"If a file name is provided, the applied forces for all actuators will be constructed "
		"from the forces_file instead of from the states.  This option should be used "
		"to calculated joint loads from static optimization results.");
	_propertySet.append(&_forcesFileNameProp);

	_jointNamesProp.setName("joint_name");
	_jointNamesProp.setComment("Name of the joint on which to perform the analysis."
		"This should be the glenohumeral joint (in Open Shoulder Model: glenohumeral.");
	_propertySet.append(&_jointNamesProp);

	_onBodyProp.setName("apply_on_body");
	_onBodyProp.setComment("Choice of body (parent or child) for which the reaction "
		"load is calculated.  Child body is default.");
	_propertySet.append(&_onBodyProp);

	_inFrameProp.setName("express_in_frame");
	_inFrameProp.setComment("Choice of frame (ground, parent, or child) in which the calculated "
		"reaction is expressed.  Parent body is default.");
	_propertySet.append(&_inFrameProp);

	_athetaProp.setName("theta_at_glenoid_rim");
	_athetaProp.setComment("Angle along the major axis of the ellipse representing the glenoid fossa.");
	_propertySet.append(&_athetaProp);

	_aphiProp.setName("phi_at_glenoid_rim");
	_aphiProp.setComment("Angle along the minor axis of the ellipse representing the glenoid fossa.");
	_propertySet.append(&_aphiProp);
}

//=============================================================================
// CONSTRUCTION METHODS
//=============================================================================
//_____________________________________________________________________________
/**
 * Setup the ReactionList which controls the joints included in the analysis,
 * the bodies the loads are applied to, and the frames the loads are expressed
 * in.
 */
void JointReaction_GH::setupReactionList()
{
	/* check length of property arrays.  if one is empty, set it to default*/
	if(_jointNames.getSize() == 0) {
		cout << "\nNo joint is specified in joint_name.  Setting to ALL.\n";
		_jointNames.setSize(1);
		_jointNames[0] = "ALL";}
	if (_onBody.getSize() == 0) {
		cout << "\nNo body is specified in apply_on_body.  Setting to parent\n";
		_onBody.setSize(1);
		_onBody[0] = "parent";
	}
	if (_inFrame.getSize() == 0) {
		cout << "\nNo body is specified in express_in_frame.  Setting to child\n";
		_inFrame.setSize(1);
		_inFrame[0] = "child";
	}

	/* get the joint set and  body set from the dynamics engine model*/
	const JointSet& jointSet = _model->getJointSet();
	const BodySet& bodySet = _model->getBodySet();
	int numJoints = jointSet.getSize();
	int numBodies = bodySet.getSize();

	/* get the ground body index in the body set*/
	int groundIndex = bodySet.getIndex("ground", 0);

	/* check if jointNames is specified to "ALL".  if yes, setup to 
	*  compute reactions loads for all joints*/
	std::string firstNameEntry = _jointNames.get(0);
	// convert to upper case
	std::transform(firstNameEntry.begin(),firstNameEntry.end(),firstNameEntry.begin(), ::toupper);
	if (firstNameEntry == "ALL") {
		_jointNames.setSize(numJoints);
		for(int i=0;i<numJoints;i++) {
			Joint& joint = jointSet.get(i);
			_jointNames.set(i, joint.getName());
		}
	}
	int numJointNames = _jointNames.getSize();

	/* check that _jointNames Array is of length 1.  If not, set lengths to default of 1 and
	*  set the values so that all reactions will be reported on the child body, expressed
	*  in the child frame.*/
	if (_onBody.getSize() == 1);
	else {
		cout << "\n WARNING: apply_on_body list is not of length 1."
			<< "\n Reaction load will be reported on the child body.\n";
		_onBody.setSize(1);
		_onBody[0] = "child";
	}

	if (_inFrame.getSize() == 1);
	else {
		cout << "\n WARNING: express_in_frame list is not of length 1."
			<< "\n Reaction load will be reported in the child body frame.\n";
		_inFrame.setSize(1);
		_inFrame[0] = "child";
	}

	/* setup the JointReaction_GHKey and, for valid joint names, determine and set the 
	*  reactionIndex, onBodyIndex, and inFrameIndex of each JointReaction_GHKey */

	_reactionList.setSize(0);
	int listNotEmptyFlag = 0;
	if (_jointNames.getSize() > 1){
		cout << "\nWARNING: The analysis requires only one joint (glenohumeral)";
		cout << "\nThe first joint on the list will be chosen for analysis.\n";
	}
		
	JointReaction_GHKey currentJoint;
	int validJointFlag = 0;
	for (int j=0; j<numJoints; j++) {
		Joint& joint = jointSet.get(j);
		if (_jointNames.get(0) == joint.getName()) {
			validJointFlag++;
			listNotEmptyFlag++;
			currentJoint.jointName = joint.getName();
			std::string childName = joint.getBody().getName();
			int childIndex = bodySet.getIndex(childName, 0);
			std::string parentName = joint.getParentBody().getName();
			int parentIndex = bodySet.getIndex(parentName, 0);

			/* set index that correponds to the appropriate index of the 
			*  computeReactions arguements forcesVec and momentsVec.*/
			currentJoint.reactionIndex = childIndex;

			/* set the onBodyIndex to the parent body*/
			currentJoint.onBodyIndex = parentIndex;
				
			/* set the inFrameIndex to the parent*/
			currentJoint.inFrameIndex = parentIndex;

			/* set the onBodyIndex to either the parent or child body*/
			std::string whichBody = _onBody[0];

			//convert whichBody to lower case
			std::transform(whichBody.begin(), whichBody.end(), whichBody.begin(), ::tolower);

			if (whichBody == "parent") {
				currentJoint.onBodyIndex = parentIndex;
			}
			else if (whichBody == "child") { currentJoint.onBodyIndex = childIndex; }
			else {
				currentJoint.onBodyIndex = childIndex;
				cout << "\nWARNING:  " << whichBody << " is not a valid choice for apply_on_body";
				cout << "\nSetting to apply " << currentJoint.jointName << " load to the child body.\n";
			}

			/* set the inFrameIndex to either the ground, child, or parent*/
			std::string whichFrame = _inFrame[0];

			// convert to lower case
			std::transform(whichFrame.begin(), whichFrame.end(), whichFrame.begin(), ::tolower);

			if (whichFrame == "child") {
				currentJoint.inFrameIndex = childIndex;
			}
			else if (whichFrame == "parent") {
				currentJoint.inFrameIndex = parentIndex;
			}
			else if (whichFrame == "ground") { currentJoint.inFrameIndex = groundIndex; }
			else {
				currentJoint.inFrameIndex = groundIndex;
				cout << "\nWARNING:  " << whichFrame << " is not a valid choice for express_in_frame";
				cout << "\nSetting to express " << currentJoint.jointName << " load in the ground frame.\n";
			}

			_reactionList.append(currentJoint);
			break;
		}
				
	}

	if (validJointFlag == 0)
		cout << "\nWARNING: " << _jointNames.get(0) << " is not a valid joint.\n";

	if(listNotEmptyFlag ==0) {
		cout << "Setting up _reactionList to include first joint in the model.\n";
		_jointNames.setSize(1);
		_jointNames[0] = jointSet.get(0).getName();
		setupReactionList();
	}
}
	



//_____________________________________________________________________________
/**
 * Construct a description for the joint reaction loads files.
 */
void JointReaction_GH::
constructDescription()
{
	string descrip;

	descrip = "\nThis file contains the GH reaction force and stability value.\n";
	descrip += "The stability value is 0 in the middle of the glenoid ";
	descrip += "and 1 at the rim of the glenoid.\n";
	descrip += "\nUnits are S.I. units (seconds, meters, Newtons, ...)";

	setDescription(descrip);
	
}

//_____________________________________________________________________________
/**
 * Construct column labels for the output results.
 *
 * For analyses that run during a simulation, the first column is almost
 * always time.  For the purpose of example, the code below adds labels
 * appropriate for recording the translation and orientation of each
 * body in the model.
 *
 * This method needs to be called as necessary to update the column labels.
 */
void JointReaction_GH::
constructColumnLabels()
{
	if(_model==NULL) return;

	Array<string> labels;
	labels.append("time");

	const BodySet& bodySet = _model->getBodySet();

	//  For the joint listed in _reactionList, append 3 column labels for forces
	std::string jointName = _reactionList.get(0).jointName;
	std::string onBodyName = bodySet.get(_reactionList.get(0).onBodyIndex).getName();
	std::string inFrameName = bodySet.get(_reactionList.get(0).inFrameIndex).getName();
	std::string labelRoot = jointName + "_on_" + onBodyName + "_in_" + inFrameName;
	labels.append(labelRoot + "_fx");
	labels.append(labelRoot + "_fy");
	labels.append(labelRoot + "_fz");

	// Add a label for the GH stability value
	labels.append("GH_stability");

	setColumnLabels(labels);
}

//_____________________________________________________________________________
/**
 * Load actuation storage from file.
 *
 * If called, this method sets _storeActuation to the
 * forces data in _forcesFileName
 */
void JointReaction_GH::
loadForcesFromFile()
{
	delete _storeActuation; _storeActuation = NULL;
	// check if the forces storage file name is valid and, if so, load the file into storage
	if(_forcesFileNameProp.isValidFileName()) {
		
		cout << "\nLoading actuator forces from file " << _forcesFileName << "." << endl;
		_storeActuation = new Storage(_forcesFileName);
		int storeSize = _storeActuation->getSmallestNumberOfStates();
		
		cout << "Found " << storeSize << " actuator forces with time stamps ranging from "
			<< _storeActuation->getFirstTime() << " to " << _storeActuation->getLastTime() << "." << endl;

		// check if actuator set and forces file have the same actuators
		bool _containsAllActuators = true;
		int actuatorSetSize = _model->getActuators().getSize();
		if(actuatorSetSize > storeSize){
			cout << "The forces file does not contain enough actuators." << endl;
			_containsAllActuators = false;
		}
		else {
			for(int actuatorIndex=0;actuatorIndex<actuatorSetSize;actuatorIndex++)
			{
				std::string actuatorName = _model->getActuators().get(actuatorIndex).getName();
				int storageIndex = _storeActuation->getStateIndex(actuatorName,0);
				if(storageIndex == -1) {
					cout << "\nThe actuator " << actuatorName << " was not found in the forces file." << endl;
					_containsAllActuators = false;
				}
			}
		}

		if(_containsAllActuators) {
			if(storeSize> actuatorSetSize) cout << "\nWARNING:  The forces file contains actuators that are not in the model's actuator set." << endl;
			_useForceStorage = true;
			cout << "WARNING:  Ignoring fiber lengths and activations from the states since " << _forcesFileNameProp.getName() << " is also set." << endl;
			cout << "Actuator forces will be constructed from " << _forcesFileName << "." << endl;
		}
		else {
			_useForceStorage = false;
			cout << "Actuator forces will be constructed from the states." << endl;
		}
	}

	else {
		cout << "WARNING:  " << _forcesFileNameProp.getName() << " is not a valid file name." << endl;
		cout << "Actuator forces will be constructed from the states." << endl;
		_useForceStorage = false;
	}
}

//_____________________________________________________________________________
/**
 * Set up storage objects.
 *
 * The storage objects in the analysis are used to record
 * the results of the analysis and write them to file.  
 */
void JointReaction_GH::
setupStorage()
{
	// Reaction Loads
	_storeReactionLoads.reset(0);
	_storeReactionLoads.setName("Glenohumeral Joint Force and Stability");
	_storeReactionLoads.setDescription(getDescription());
	_storeReactionLoads.setColumnLabels(getColumnLabels());

	// Actuator forces - if a forces file is specified, load the forces storage data to _storeActuation
	if(!(_forcesFileName == "")) loadForcesFromFile();

}


//=============================================================================
// GET AND SET
//=============================================================================
//_____________________________________________________________________________
/**
 * Set the model for which this analysis is to be run.
 *
 * Sometimes the model on which an analysis should be run is not available
 * at the time an analysis is created.  Or, you might want to change the
 * model.  This method is used to set the model on which the analysis is
 * to be run.
 *
 * @param aModel Model pointer
 */
void JointReaction_GH::
setModel(Model& aModel)
{
	// SET THE MODEL IN THE BASE CLASS
	Analysis::setModel(aModel);

	// UPDATE VARIABLES IN THIS CLASS
	setupReactionList();
	constructDescription();
	constructColumnLabels();
	//setupStorage();
	_dydt.setSize(_model->getNumStateVariables());
	int numJoints = _reactionList.getSize();
	// set size of working array of loads. 
	// 3 for force, and 1 for the GH stability value
	_Loads.setSize(4);
}


//=============================================================================
// ANALYSIS
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute and record the results.
 *
 * This method computes the reaction loads at all joints in the model, then
 * truncates the results to contain only the loads at the requested joints
 * and finally, if necessary, modifies the loads to be acting on the specified
 * body and expressed in the specified frame
 *
 * @param aT Current time in the simulation.
 * @param aX Current values of the controls.
 * @param aY Current values of the states.
 */
int JointReaction_GH::
record(const SimTK::State& s)
{
	/** if a forces file is specified replace the computed actuation with the 
	    forces from storage.*/
	SimTK::State s_analysis = s;

	_model->updMultibodySystem().realize(s_analysis, s.getSystemStage());
	if(_useForceStorage){

		const Set<Actuator> *actuatorSet = &_model->getActuators();
		int nA = actuatorSet->getSize();
		Array<double> forces(0,nA);
		_storeActuation->getDataAtTime(s.getTime(),nA,forces);
		int storageIndex = -1;
		for(int actuatorIndex=0;actuatorIndex<nA;actuatorIndex++)
		{
			//Actuator* act = dynamic_cast<Actuator*>(&_forceSet->get(actuatorIndex));
			std::string actuatorName = actuatorSet->get(actuatorIndex).getName();
			storageIndex = _storeActuation->getStateIndex(actuatorName, 0);
			if(storageIndex == -1){
				cout << "The actuator, " << actuatorName << ", was not found in the forces file." << endl;
				break;
			}
			actuatorSet->get(actuatorIndex).overrideForce(s_analysis,true);
			actuatorSet->get(actuatorIndex).setOverrideForce(s_analysis,forces[storageIndex]);
		}
	}

	// VARIABLES
	int numBodies = _model->getNumBodies();

	/** define 2 variable length vectors of Vec3 vectors to contain calculated  
	*   forces and moments for all the bodies in the model */
	Vector_<Vec3> allForcesVec(numBodies);
	Vector_<Vec3> allMomentsVec(numBodies);

	//// BodySet and JointSet and ground body index
	const BodySet& bodySet = _model->getBodySet();
	const JointSet& jointSet = _model->getJointSet();
	Body &ground = _model->getSimbodyEngine().getGroundBody();

	/* Calculate All joint reaction forces and moments.
	*  Applied to child bodies, expressed in ground frame.  
	*  computeReactions realizes to the acceleration stage internally
	*  so you don't have to call realize in this analysis.*/ 
	_model->getSimbodyEngine().computeReactions(s_analysis, allForcesVec, allMomentsVec);

	/* retrieve desired joint reaction, and convert
	*  to parent reference frame (scapula)*/
	JointReaction_GHKey currentKey = _reactionList[0];
	const Joint& joint = jointSet.get(currentKey.jointName);
	Vec3 force = -allForcesVec[currentKey.reactionIndex];

	/* Take glenoid orientation into account */
	double Rgt_x[3] = { 0.9511, -0.2724, -0.1455 };
	double Rgt_y[3] = { 0.2693, 0.9622, -0.0412 };
	double Rgt_z[3] = { 0.1513, 0, 0.9885 };
//	double Rgt_x[3] = { 1, 0, 0 };
//	double Rgt_y[3] = { 0, 1, 0 };
//	double Rgt_z[3] = { 0, 0, 1 };

	double gforce[3] = { 0 };
	/* The glenoid orientation matrix came from the original Delft model
	* which had x pointing laterally, y superiorly and z posteriorly */
	gforce[0] = Rgt_x[0] * force[2] + Rgt_x[1] * force[1] - Rgt_x[2] * force[0];
	gforce[1] = Rgt_y[0] * force[2] + Rgt_y[1] * force[1] - Rgt_y[2] * force[0];
	gforce[2] = Rgt_z[0] * force[2] + Rgt_z[1] * force[1] - Rgt_z[2] * force[0];

	double phi, theta;
	double GHstab = 0;

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

	GHstab = pow((theta / _atheta), 2) + pow((phi / _aphi), 2);

	if (currentKey.onBodyIndex == currentKey.reactionIndex)
		/*Take reaction load from parent and apply on child*/
		force = -force;

	/* express load in the desired reference frame*/
	Body& expressedInBody = bodySet.get(currentKey.inFrameIndex);
	_model->getSimbodyEngine().transform(s_analysis, ground, force, expressedInBody, force);

	/* fill out row construction array*/
	for(int j=0;j<3;j++) {
		_Loads[j] = force[j];
	}
	// stability value
	_Loads[3] = GHstab;

	/* Write the reaction data to storage*/
	_storeReactionLoads.append(s.getTime(),_Loads.getSize(),&_Loads[0]);


	return(0);
}
//_____________________________________________________________________________
/**
 * This method is called at the beginning of an analysis so that any
 * necessary initializations may be performed.
 *
 * This method is meant to be called at the begining of an integration 
 *
 * @param s reference to the current state
 *
 * @return -1 on error, 0 otherwise.
 */
int JointReaction_GH::
begin(SimTK::State& s)
{
	if(!proceed()) return(0);
	// Read forces file here rather than during initialization
	setupStorage();

	// RESET STORAGE
	_storeReactionLoads.reset(s.getTime());

	// RECORD
	int status = 0;
	if(_storeReactionLoads.getSize()<=0) {
		status = record(s);
	}

	return(status);
}
//_____________________________________________________________________________
/**
 * This method is called to perform the analysis.  It can be called during
 * the execution of a forward integrations or after the integration by
 * feeding it the necessary data.
 *
 *
 * @param s reference to the current stateaClientData General use pointer for sending in client data.
 *
 * @return -1 on error, 0 otherwise.
 */
int JointReaction_GH::
step( const SimTK::State& s, int stepNumber)
{
	if(!proceed(stepNumber)) return(0);

	record(s);

	return(0);
}
//_____________________________________________________________________________
/**
 * This method is called at the end of an analysis so that any
 * necessary finalizations may be performed.
 *
 * @param s reference to the current state
 *
 * @return -1 on error, 0 otherwise.
 */
int JointReaction_GH::
end(SimTK::State& s)
{
	if(!proceed()) return(0);

	record(s);

	return(0);
}




//=============================================================================
// IO
//=============================================================================
//_____________________________________________________________________________
/**
 * Print results.
 * 
 * The file names are constructed as
 * aDir + "/" + aBaseName + "_" + ComponentName + aExtension
 *
 * @param aDir Directory in which the results reside.
 * @param aBaseName Base file name.
 * @param aDT Desired time interval between adjacent storage vectors.  Linear
 * interpolation is used to print the data out at the desired interval.
 * @param aExtension File extension.
 *
 * @return 0 on success, -1 on error.
 */
int JointReaction_GH::
printResults(const string &aBaseName,const string &aDir,double aDT,
				 const string &aExtension)
{
	// Reaction Loads
	Storage::printResult(&_storeReactionLoads,aBaseName+"_"+getName()+"_ReactionLoads",aDir,aDT,aExtension);

	return(0);
}


