#ifndef actuatorSurfaceSource_H
#define actuatorSurfaceSource_H

#include "List.H"
#include "PtrList.H"
#include "dictionary.H"
#include "floatScalar.H"
#include "fvMatrix.H"
#include "fvMesh.H"
#include "scalar.H"
#include "slicedVolFieldsFwd.H"
#include "typeInfo.H"
#include "vector.H"
#include "actuatorSurfaceElement.H"
#include "cellSetOption.H"
#include "volFieldsFwd.H"


namespace Foam
{
namespace fv 
{

class actuatorSurfaceSource 
:
  public cellSetOption
{
  protected:
    
    //protected data
      
      // total number of points that make up the actuatorSurfaceSource
      label nPoints_;

      vector freeStreamVelocity_;

      vector freeStreamDirection_;

      // total force field due to the actuator surface
      volVectorField forceField_;
      
      // surface comprising of Surfaceelements 
      PtrList<actuatorSurfaceSource> elements_;
    
    // protected member functions
      
      // create the surface that is going to make up the surface of the nacelle 
      void createElements();

      //- read dictionary
      bool read(const dictionary &dict);

      // create output fil
      virtual void createOutputFile();

      // write relevant metrics needed for a nacelle flow
      void writePerf();
  public:
    // public methods
      
      TypeName("actuatorSurfaceSource");

      // return a reference to the selected fvOption Model
      static autoPtr<actuatorSurfaceSource> New
      (
       const word &name,
       const dictionary &dict,
       const fvMesh &mesh
      );

      // destructor
      virtual ~actuatorSurfaceSource();
      
      //methods
        
        //access
        
          // return const reference to  total force vector
          const vector &force();
          
          // return const reference to the force field
          const volVectorField &forceField();
          
          //  return reference to the element pointer list that makes up the AS
          PtrList<actuatorSurfaceSource> &elements();

        //edit

          // rotate the surface about a specified axis
          void rotate(vector rotationPoint, vector axis, scalar radians);
          
          void setSpeed(vector point, vector axis, scalar radians);


        // source term addition
          
          // source term to momentum equation
          virtual void addSup
          (
           fvMatrix<vector>&eqn,
           const label fieldI
          );
          
          // source terms to turbulence scalar transport
          virtual void addSup
          (
           fvMatrix<scalar>& eqn,
           const label fieldI
          );
        








};

} // end of namespace fv 
} // end of namespace Foam

#endif // !actuatorSurfaceSource_H
