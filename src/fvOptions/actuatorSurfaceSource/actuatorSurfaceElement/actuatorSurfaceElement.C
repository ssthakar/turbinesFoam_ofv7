#include "actuatorSurfaceElement.H"
#include "dictionary.H"
#include "fvMesh.H"

// * * * * *  * * * * * * * Protected Methods * * * * * * * * * * * * * * * //

// 4 point cosine function 
scalar Foam::fv::actuatorSurfaceElement::phiR(scalar r)
{
    scalar phiR = 0.0;
    scalar pi_ = Foam::constant::mathematical::pi;
    if(r<=1.5)
    {
        phiR = 0.25 + (sin(0.25*pi_*(2*r+1)) - sin(0.25*pi_*(2*r+1)))/(2*pi_);
    }
    else if(r > 1.5 && r < 2.5)
    {
        phiR = 5.0/8 - (sin(0.25*pi_*(2*r+1))/(2*pi_));
    }
    else
    {
        phiR = 0.0;
    }
    return phiR;
}

// discrete delta function to calculate flow velocity at element location
scalar Foam::fv::actuatorSurfaceElement::discreteDeltaFunction(label cellI)
{ 
    scalar delta_ = 0.0;
    scalar mesh_length_scale = Foam::cbrt(mesh_.V()[cellI]); // this is heavily based on the assumption that the cells in the cell set are cubes 
    vector r = mesh_.C()[cellI] - position_;
    delta_ = phiR(r.x()/mesh_length_scale)*phiR(r.y()/mesh_length_scale)*phiR(r.z()/mesh_length_scale);
    delta_ = delta_/mesh_.V()[cellI];
    return delta_;
}

// discrete delta function to calculate the flow velocity at any given location
scalar Foam::fv::actuatorSurfaceElement::discreteDeltaFunction(label cellI,const point position)
{
    scalar delta = 0.0; // intialize the delta value 
    scalar mesh_length_scale = Foam::cbrt(mesh_.V()[cellI]);
    vector r = mesh_.C()[cellI] - position;
    delta = phiR(r.x()/mesh_length_scale)*phiR(r.y()/mesh_length_scale)*phiR(r.z()/mesh_length_scale);
    delta = delta/mesh_.V()[cellI];
    return delta;
}


// calculate the inflow velocity at element position from background fvMesh
void Foam::fv::actuatorSurfaceElement::calculateInflowVelocity(const volVectorField &Uin)
{
   inflowVelocity_ = vector(0.0,0.0,0.0);
   forAll(mesh_.cells(),cellI)
   {
    inflowVelocity_ += Uin[cellI]*discreteDeltaFunction(cellI)*mesh_.V()[cellI];
   } 
}

vector Foam::fv::actuatorSurfaceElement::calculateInflowVelocity(const volVectorField &Uin, vector inflowPosition)
{
    vector inflowVelocity = vector(0.0,0.0,0.0); // intialize flow velocity as 0.
    label posCellI = findCell(inflowPosition); // get the cell that has the position
    forAll(mesh_.cells(),cellI)
    {
        inflowVelocity += Uin[cellI]*discreteDeltaFunction(cellI,inflowPosition)*mesh_.V()[cellI];
    }
    return inflowVelocity;
}

// function to find the cell that holds the point
Foam::label Foam::fv::actuatorSurfaceElement::findCell
(
    const point& location
)
{
    if (Pstream::parRun()) //check if code is running in parallel
    {
        if (meshBoundBox_.containsInside(location)) //check for point in processor bounding box
        {
            if (debug) 
            {
                Pout<< "Looking for cell containing " << location
                    << " inside bounding box:" << endl
                    << meshBoundBox_ << endl;
            }
            return mesh_.findCell(location); // return cell Label
        }
        else
        {
            if (debug)
            {
                Pout<< "Cell not inside " << meshBoundBox_ << endl;
            }
            return -1;
        }
    }
    else
    {
        return mesh_.findCell(location);;
    }
}

// calculate normal force acting on actuator surface 
void Foam::fv::actuatorSurfaceElement::calculateForce(const volVectorField &Uin)
{
    // calculate inflow velocity at actuator surface element centre 
    calculateInflowVelocity(Uin);
    // current time stepientu
    scalar deltaT = time_.deltaT().value();

    // calculate the normal force acting on the actuator surface element from direct immersed boundary method
    forceNormal_ = lengthScale_*((inflowVelocity_ - velocity_)&planformNormal_)*planformNormal_/(deltaT);

    // calculate the tangential force acting on the actuator surface element from turbulennt coeffic
    forceTangential_ = 0.5*turbulentFrictionCoeff(tangentialVelocity_,lengthScale_)*mag(freeStreamVelocity_)*mag(freeStreamVelocity_)*tangentialDir_;
}


// needs reference to viscosity of the current simulation to calculate the Reynolds number
scalar Foam::fv::actuatorSurfaceElement::turbulentFrictionCoeff(vector &velocityAtPoint, scalar &UpStreamLengthScale)
{
    // Reynolds number based on the length scale and inflow velocity sampled at sampling point
    scalar Re_x=0.0;
    Re_x = mag(inflowVelocity_)*UpStreamLengthScale/nu_;
    return 0.37*pow(log(Re_x),-2.584);
} 

//  Constructors

Foam::fv::actuatorSurfaceElement::actuatorSurfaceElement
(
    const dictionary &dict,
    const fvMesh &mesh,
    const Time &time
)
:
dict_(dict),
mesh_(mesh),
time_(time),
meshBoundBox_(mesh_.points(),false) // intialize the meshBoundBox from pointField, false is to not reduce over processors

{

}























