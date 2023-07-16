#ifndef actuatorSurfaceElement_H
#define actuatorSurfaceElement_H

#include "OFstream.H"
#include "autoPtr.H"
#include "boundBox.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicStallModel.H"
#include "interpolationCellPoint.H"
#include "profileData.H"
#include "addedMassModel.H"

namespace Foam
{
namespace fv 
{

class actuatorSurfaceElement 
{
protected:
  
  //protected data 
    
    //- because the class members are references, they need to be intialized when the constructor is called
    //- this can be done through an intializer list

      //- dictionary to store in data related to the actuator surface 
      const dictionary dict_;

      //- reference to the mesh object
      const fvMesh &mesh_;

      // Time object reference 
      const Time &time_; 

      //- Mesh bounding box 
      boundBox meshBoundBox_;

      //- location of the element, used in source 
      vector position_;

      //- Cell Id of the cell containing the actuator surface element
      label PosCellI_;

      //- sampling position for calculating Vrel  
      vector samplingPosition_;

      //- position to use to calculate the tangential sampling point
      vector tangentialVelSamplingPoint_;

      //- tangential Velocity at tangential Velocity Sampling point
      vector tangentialVelocity_;

      //- direction of the tangential force that is acting on the surface element
      vector tangentialDir_;

      //- lagrangian velocity of the element 
      vector velocity_;

      //- Eulerian Free stream velocity
      vector freeStreamVelocity_;
    
      //- FreeStream Direction
      vector freeStreamDirection_;

      //- the normal vector to the surface
      vector planformNormal_;

      //- Relative velocity 
      vector relativeVelocity_;

      //- inflow velocity (Eulerian)
      vector inflowVelocity_;
    
      //- Fluid force on element along normal to surface
      vector forceTangential_;

      //- Fluid force along the tangential direction
      vector forceNormal_;

      //- angle between normal vector to surface and inflow Eulerian 
      //- eulerian velocity
      scalar angleOfAttack_;

      //- don't know what the fuck this is or why the fuck this is
      scalar angleOfAttackGeom_;

      //- Reynolds number based on the inflow velocity and friction lenght
      scalar Re_;

      // - angular velocity of the turbine in rpm
      scalar omega_;

      //- viscosity from transport properties
      scalar nu_;

      //- length scale based on the mesh
      scalar lengthScale_;

      //-  switch to write output data
      bool writePerf_;
      
      word samplingMethod_;  

      //- output stream pointer
      OFstream *outputFile_;

  //Protected Member Functions
    //- Rotate a vector
      void rotateVector
      (
      vector &vectorToRotate,
      vector rotationPoint,
      vector axis,
      scalar radians
      );

      //- Find cell containing location
      label findCell(const point &location);

      // four point smooth cosine function
      scalar phiR(scalar r);

      //- discreteDeltaFunction to sample and source stuff
      scalar discreteDeltaFunction(label cellI);

      //- discrete deltaFunction for sampling velocity at a particular location
      scalar discreteDeltaFunction(label cellI,const point samplingPosition);

      //
      void multiplyForceRho(const volScalarField &rho);

      //
      void applyForceField(volVectorField &forceField);

      //-calculate velocity at element location
      void calculateInflowVelocity(const volVectorField &Uin);
      
      //- calcualte the velocity given position
      vector calculateInflowVelocity(const volVectorField & Uin, vector inflowPosition);

      scalar turbulentFrictionCoeff(vector &velocityAtPoint,scalar &lengthScale);

      // set the lagrangian velocity of the actuator surface element based 
      // on the TSR of the turbine

      
public:

    // - constructor 
    actuatorSurfaceElement
    (
     const dictionary &dict,
     const fvMesh &mesh,
     const Time &time
    );
    
    // set element speed based on turbine TSR 
    void velocityFromTSR();

    // calculate the normal and tangential Force acting on the surfaceElement
    void calculateForce(const volVectorField &Uin);
    
};


} //end of namespace  fv 
} //end of namespace Foam


#endif // !actuatorsSurfaceElement_H
