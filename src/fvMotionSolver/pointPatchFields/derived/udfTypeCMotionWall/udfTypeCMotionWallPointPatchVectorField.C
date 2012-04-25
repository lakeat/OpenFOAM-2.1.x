/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "udfTypeCMotionWallPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar udfTypeCMotionWallPointPatchVectorField::pointGlobalAngle
(
    const point& globalOrigin,
    const point& globalPoint
) const
{
    // Assume the global Y axis unit vector is (0,1,0)
    vector tmpVector = globalPoint-globalOrigin;
    scalar theta = 0;
    if ( mag(tmpVector.x()) > VSMALL )
    {
        scalar tmpTheta = radToDeg(atan(mag(tmpVector.y())/mag(tmpVector.x())));
        if (tmpVector.x()>=0 && tmpVector.y()>=0)
        {
            theta = tmpTheta;
        }
        else if (tmpVector.x()<0 && tmpVector.y()>=0)
        {
            theta = 180-tmpTheta;
        }
        else if (tmpVector.x()<=0 && tmpVector.y()<0)
        {
            theta = tmpTheta+180;
        }
        else
        {
            theta = 360-tmpTheta;
        }
    }
    else
    {
        if ( tmpVector.y() > 0 )
        {
            theta = 90;
        }
        else
        {
            theta = 270;
        }
    }
    
    return theta;
}


const point udfTypeCMotionWallPointPatchVectorField::twoDGloToLocTransSpecial
(
    const point& globalOrigin,
    const vector& globalYUnit,
    const point& lockEndLt,
    const point& lockEndRt,
    const point& globalPoint
) const
{
    vector tmpVector = lockEndRt-lockEndLt;
    point lockEndCt = 0.5*(lockEndRt+lockEndLt);
    tensor rotateTmp1(0,-1,0,1,0,0,0,0,1);
    vector localYUnit = rotateTmp1 & (tmpVector/mag(tmpVector));
    vector translationVector = globalOrigin-lockEndCt;
    
    vector n1=localYUnit;
    vector n2=globalYUnit/mag(globalYUnit);
    scalar s = n1 & n2;
    vector n3 = n1 ^ n2;
    tensor rotationMatrix = s*I
                          + (1 - s)*sqr(n3)/(magSqr(n3) + VSMALL)
                          + (n2*n1 - n1*n2);
    
    point temp = rotationMatrix & (globalPoint + translationVector);
    temp.z() = globalPoint.z();
    
    return temp;
}


const point udfTypeCMotionWallPointPatchVectorField::twoDLocToGloTransSpecial
(
    const point& globalOrigin,
    const vector& globalYUnit,
    const point& lockEndLt,
    const point& lockEndRt,
    const point& localPoint
) const
{
    vector tmpVector = lockEndRt-lockEndLt;
    point lockEndCt = 0.5*(lockEndRt+lockEndLt);
    tensor rotateTmp1(0,-1,0,1,0,0,0,0,1);
    vector localYUnit = rotateTmp1 & (tmpVector/mag(tmpVector));
    vector translationVector = lockEndCt-globalOrigin;
    
    vector n1=globalYUnit/mag(globalYUnit);
    vector n2=localYUnit;
    scalar s = n1 & n2;
    vector n3 = n1 ^ n2;
    tensor rotationMatrix = s*I
                          + (1 - s)*sqr(n3)/(magSqr(n3) + VSMALL)
                          + (n2*n1 - n1*n2);
    
    point temp = (rotationMatrix & localPoint) + translationVector;
    temp.z() = localPoint.z();
    
    return temp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

udfTypeCMotionWallPointPatchVectorField::
udfTypeCMotionWallPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    initialPoints_(p.localPoints()),
    lockEndPoints_(),
    ctrlPoints_(),
    staChangingTime_(0.0),
    endChangingTime_(0.0),
    globalOrigin_(vector::zero),
    globalYUnit_(vector::zero),
    endBCType_(0.0),
    endBCTypeValL_(0.0),
    endBCTypeValR_(0.0)
{}


udfTypeCMotionWallPointPatchVectorField::
udfTypeCMotionWallPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    staChangingTime_(readScalar(dict.lookup("staChangingTime"))),
    endChangingTime_(readScalar(dict.lookup("endChangingTime"))),
    globalOrigin_(dict.lookup("globalOrigin")),
    globalYUnit_(dict.lookup("globalYUnit")),
    endBCType_(readScalar(dict.lookup("endBCType"))),
    endBCTypeValL_(readScalar(dict.lookup("endBCTypeValL"))),
    endBCTypeValR_(readScalar(dict.lookup("endBCTypeValR")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("lockEndPoints"))
    {
        lockEndPoints_ = vectorField("lockEndPoints", dict , 8);
    }

    if (dict.found("ctrlPoints"))
    {
        ctrlPoints_ = vectorField("ctrlPoints", dict , 17);
    }

    if (dict.found("initialPoints"))
    {
        initialPoints_ = vectorField("initialPoints", dict , p.size());
    }
    else
    {
        initialPoints_ = p.localPoints();
    }
}


udfTypeCMotionWallPointPatchVectorField::
udfTypeCMotionWallPointPatchVectorField
(
    const udfTypeCMotionWallPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    initialPoints_(ptf.initialPoints_, mapper),
    lockEndPoints_(ptf.lockEndPoints_),
    ctrlPoints_(ptf.ctrlPoints_),
    staChangingTime_(ptf.staChangingTime_),
    endChangingTime_(ptf.endChangingTime_),
    globalOrigin_(ptf.globalOrigin_),
    globalYUnit_(ptf.globalYUnit_),
    endBCType_(ptf.endBCType_),
    endBCTypeValL_(ptf.endBCTypeValL_),
    endBCTypeValR_(ptf.endBCTypeValR_)
{}


udfTypeCMotionWallPointPatchVectorField::
udfTypeCMotionWallPointPatchVectorField
(
    const udfTypeCMotionWallPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    initialPoints_(ptf.initialPoints_),
    lockEndPoints_(ptf.lockEndPoints_),
    ctrlPoints_(ptf.ctrlPoints_),
    staChangingTime_(ptf.staChangingTime_),
    endChangingTime_(ptf.endChangingTime_),
    globalOrigin_(ptf.globalOrigin_),
    globalYUnit_(ptf.globalYUnit_),
    endBCType_(ptf.endBCType_),
    endBCTypeValL_(ptf.endBCTypeValL_),
    endBCTypeValR_(ptf.endBCTypeValR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void udfTypeCMotionWallPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchVectorField::autoMap(m);

    initialPoints_.autoMap(m);
}


void udfTypeCMotionWallPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const udfTypeCMotionWallPointPatchVectorField& uSDoFptf =
    refCast
    <
        const udfTypeCMotionWallPointPatchVectorField
    >(ptf);

    fixedValuePointPatchVectorField::rmap(uSDoFptf, addr);

    initialPoints_.rmap(uSDoFptf.initialPoints_, addr);
}


void udfTypeCMotionWallPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& ptPatch = this->patch();
    
//------------------------------------------------------------------------------
// Step (1) Get the bc points and their locations before deformation
//------------------------------------------------------------------------------

    Foam::Time startingTime
    (
        Foam::Time::controlDictName,
        t.rootPath(),
        t.caseName()
    );
    startingTime.setTime(0,0);
    Foam::fvMesh meshOrg
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            startingTime.timeName(),
            startingTime,
            Foam::IOobject::MUST_READ
        )
    );
    const pointField& pointsOrg = meshOrg.points();
    label nMarkedPointsOrg = 0;
    labelList bcPtsIDsOrg;
    vectorField bcPointsOrg;
    label patchIndexOrg = meshOrg.boundaryMesh().findPatchID(ptPatch.name());
    const labelList& mpOrg = meshOrg.boundaryMesh()[patchIndexOrg].meshPoints();
    bcPtsIDsOrg.setSize(mpOrg.size());
    bcPointsOrg.setSize(mpOrg.size());
    //Info << "Moving " << mpOrg.size() << " boundary points" << endl;
    for
    (
        label pickedPoint = 0;
        pickedPoint < mpOrg.size();
        pickedPoint += 1
    )
    {
        bcPtsIDsOrg[nMarkedPointsOrg] = mpOrg[pickedPoint];
        nMarkedPointsOrg++;
    }
    
    forAll (bcPtsIDsOrg, i)
    {
        bcPointsOrg[i] = pointsOrg[bcPtsIDsOrg[i]];
    }

    //Info << "starting time : " << startingTime.value() 
    //     << "    ------->    " << tab
    //     << "current  time : " << timeNew
    //     << endl;
    
//------------------------------------------------------------------------------
// Step (2) Get the bc points and their locations (current)
//------------------------------------------------------------------------------

    const pointField& points = mesh.points();
    label nMarkedPoints = 0;
    labelList bcPtsIDsCur;
    vectorField bcPointsCur;
    label patchIndex = mesh.boundaryMesh().findPatchID(ptPatch.name());
    const labelList& mp = mesh.boundaryMesh()[patchIndex].meshPoints();
    bcPtsIDsCur.setSize(mp.size());
    bcPointsCur.setSize(mp.size());
    //Info<< "There is " << mp.size()
    //    << " grid points IN TOTAL on the surface." << endl;

    for
    (
        label pickedPoint = 0;
        pickedPoint < mp.size();
        pickedPoint += 1
    )
    {
        bcPtsIDsCur[nMarkedPoints] = mp[pickedPoint];
        nMarkedPoints++;
    }
    
    forAll (bcPtsIDsCur, i)
    {
        bcPointsCur[i] = points[bcPtsIDsCur[i]];
    }
    
    pointField bcPointsNew=bcPointsOrg;
    
//------------------------------------------------------------------------------
// Step (3) Preparation before working on the corners
//------------------------------------------------------------------------------

    //-The following parameters are defined to make life easy.
    //scalar pi = mathematicalConstant::pi;
    
    //-Whether a point is a corner points is judged by checking its angle, see
    // if it lies in between the two lock ends.
    scalarField lockEndAngles=mag(lockEndPoints_);
    for (label i=0; i<lockEndPoints_.size(); i++)
    {
        lockEndAngles[i]=pointGlobalAngle(globalOrigin_,lockEndPoints_[i]);
    }

//------------------------------------------------------------------------------
// Step (4) Working on the 1st corner
//------------------------------------------------------------------------------

//Info<< "startingTime= " << startingTime.value() << endl;
//Info<< "t= " << t.value() << endl;
//Info<< "deltaTValue= " << t.deltaTValue() << endl;
scalar deformStaTime = startingTime.value()+staChangingTime_*t.deltaTValue();
scalar deformEndTime = startingTime.value()+endChangingTime_*t.deltaTValue();
if
(
    t.value() >= deformStaTime &&
    t.value() < deformEndTime
)
{

    Info<< "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
        << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
        << endl;
    label mm=0;
    label nn=1;
    
    //-Increase the ctrl points by two, adding the two lock ends. These are the
    // real control points used by ALGLIB.
    pointField ctrlPointsTmp;
    ctrlPointsTmp.setSize(ctrlPoints_.size()+2);
    ctrlPointsTmp[0].component(vector::X)=
        -0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[0].component(vector::Y)=0;
    ctrlPointsTmp[0].component(vector::Z)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::X)=
        0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Y)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Z)=0;
    for (label i=1; i<=ctrlPoints_.size(); i++)
    {
        ctrlPointsTmp[i].component(vector::X)=
            ctrlPoints_[i-1].component(vector::X);
        ctrlPointsTmp[i].component(vector::Y)=
            ctrlPoints_[i-1].component(vector::Y);
        ctrlPointsTmp[i].component(vector::Z)=
            ctrlPoints_[i-1].component(vector::Z);
    }
    

    //-For each corner, find the corner points and then transform them to local
    // coordinates.
    //-The transformation function should be able to handle all four quadrants.
    label NCornerPoints=0;
    /*
    Info<< "lockEndAngles[0]=" << lockEndAngles[mm] << tab << tab
        << "lockEndAngles[1]=" << lockEndAngles[nn] << endl;
    */
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            //Info<< "points[" << tab << i << "]= " << bcPointsNew[i] << endl;
            NCornerPoints=NCornerPoints+1;
        }
    }
    //Info<< "There are " << NCornerPoints 
    //    << " points in ONE corner ONE slice of mesh." << endl;
    pointField locBcPointsOrg;
    locBcPointsOrg.setSize(NCornerPoints);
    pointField locBcPointsNew=locBcPointsOrg;
    scalarField delta=ctrlPointsTmp.component(vector::Y);
    
    labelList marker(bcPtsIDsCur.size(),0);
    labelList markerA(bcPtsIDsCur.size(),0);
    labelList markerB(bcPtsIDsCur.size(),0);
    labelList markerC(bcPtsIDsCur.size(),0);
    labelList markerD(bcPtsIDsCur.size(),0);
    labelList markerE(bcPtsIDsCur.size(),0);
    labelList markerF(bcPtsIDsCur.size(),0);
    labelList markerG(bcPtsIDsCur.size(),0);
    labelList markerInverse(bcPtsIDsCur.size(),0);
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        scalar ta=0;
        scalar tb=0;
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            //-Take care of the index mismatch.
            //-(1) Find the index patterns between local and global points.
            //-(2) Find the index patterns between four quadrants.
            //-(3) Find the index patterns between different slices.
            // By doing so, the assumptions should be that angle of wind attack
            // should be zero and also this should be a rectangular shape?!
            // TO-DO--------------------
            marker[i] = NCornerPoints;
            markerInverse[NCornerPoints]=i;
            NCornerPoints=NCornerPoints+1;
            ta=bcPointsNew[i].component(vector::X);
            tb=bcPointsNew[i].component(vector::Y);
        }
        forAll (bcPtsIDsCur, j)
        {
            if
            (
                (mag(bcPointsNew[j].component(vector::X)) - ta) < 1E-5 &&
                (mag(bcPointsNew[j].component(vector::Y)) - tb) < 1E-5 &&
                bcPointsNew[j].component(vector::Z) == 0
            )
            {
                if
                (
                    bcPointsNew[j].component(vector::X) < 0 &&
                    bcPointsNew[j].component(vector::Y) > 0
                )
                {
                    markerA[i]=j;
                }
                if
                (
                    bcPointsNew[j].component(vector::X) < 0 &&
                    bcPointsNew[j].component(vector::Y) < 0
                )
                {
                    markerB[i]=j;
                }
                if
                (
                    bcPointsNew[j].component(vector::X) > 0 &&
                    bcPointsNew[j].component(vector::Y) < 0
                )
                {
                    markerC[i]=j;
                };
            }
            if
            (
                (mag(bcPointsNew[j].component(vector::X)) - ta) < 1E-5 &&
                (mag(bcPointsNew[j].component(vector::Y)) - tb) < 1E-5 &&
                bcPointsNew[j].component(vector::Z) != 0
            )
            {
                if
                (
                    bcPointsNew[j].component(vector::X) > 0 &&
                    bcPointsNew[j].component(vector::Y) > 0
                )
                {
                    markerD[i]=j;
                }
                if
                (
                    bcPointsNew[j].component(vector::X) < 0 &&
                    bcPointsNew[j].component(vector::Y) > 0
                )
                {
                    markerE[i]=j;
                }
                if
                (
                    bcPointsNew[j].component(vector::X) < 0 &&
                    bcPointsNew[j].component(vector::Y) < 0
                )
                {
                    markerF[i]=j;
                }
                if
                (
                    bcPointsNew[j].component(vector::X) > 0 &&
                    bcPointsNew[j].component(vector::Y) < 0
                )
                {
                    markerG[i]=j;
                };
            }
        }
    }
    
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {            
            locBcPointsOrg[NCornerPoints]=twoDGloToLocTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                bcPointsNew[i]
            );
            /*
            Info<< "Global point " << tab << i << tab << "= "
                << "("
                << bcPointsNew[i].component(vector::X)
                << ", "
                << bcPointsNew[i].component(vector::Y)
                << ", "
                << bcPointsNew[i].component(vector::Z)
                << ")"
                << tab
                << tab
                << "Angle= " << pointGlobalAngle(globalOrigin_,bcPointsNew[i])
                << tab
                << tab
                << "Local point " << tab << NCornerPoints << tab << "= "
                << "("
                << locBcPointsOrg[NCornerPoints].component(vector::X)
                << ", "
                << locBcPointsOrg[NCornerPoints].component(vector::Y)
                << ", "
                << locBcPointsOrg[NCornerPoints].component(vector::Z)
                << ")"
                << endl;
            */
            NCornerPoints=NCornerPoints+1;
        }
    }
    
    
    //-Call ALGLIB to do the job!
    //-The job is to fit the control points with a cubic spline.
    //-This is newly added in src/, I implemented it into OpenFOAM.
    //-ALGLIB is another opensource, whoes copy right is reserved and protected
    // under GNU. DW, Nov/28/2011
    //-----------------------------
    
    //-Local x of control points are easy, they are input by the user, but local
    // y is not, we just know the delta values input by the user, so in order to
    // use the ALGLIB, we need to find the local y value and then plus the delta
    // values, so to get the final y values for the control points.
    //-Find control points' local y value, this is done by cubic spline 
    // interpolation between the local grid points.
    
    alglib::real_1d_array inputXTmp;
    alglib::real_1d_array inputYTmp;
    alglib::real_1d_array outputXTmp;
    alglib::real_1d_array outputYTmp;
    
    inputXTmp.setlength(NCornerPoints);
    inputYTmp.setlength(NCornerPoints);
    outputXTmp.setlength(ctrlPointsTmp.size());
    outputYTmp.setlength(ctrlPointsTmp.size());
    
    for(label i=0; i<NCornerPoints; i++)
    {
        inputXTmp[i]=locBcPointsOrg[i].component(vector::X);
        inputYTmp[i]=locBcPointsOrg[i].component(vector::Y);
    }
    //Info<< "inputXTmp=" << locBcPointsOrg.component(vector::X) << endl;
    //Info<< "inputYTmp=" << locBcPointsOrg.component(vector::Y) << endl;
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        outputXTmp[i]=ctrlPointsTmp[i].component(vector::X);
    }
    //Info<< "outputXTmp=" << ctrlPointsTmp.component(vector::X) << endl;
    
    alglib::ae_int_t orgSizeTmp=inputXTmp.length();
    alglib::ae_int_t newSizeTmp=outputXTmp.length();
    alglib::ae_int_t boundType = endBCType_;
    double boundTypeLeftValue = endBCTypeValL_;
    double boundTypeRightValue = endBCTypeValR_;
        
    alglib::spline1dconvcubic
    (
        inputXTmp,
        inputYTmp,
        orgSizeTmp,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputXTmp,
        newSizeTmp,
        outputYTmp
    );
    
    alglib::real_1d_array inputX;
    alglib::real_1d_array inputY;
    alglib::real_1d_array outputX;
    alglib::real_1d_array outputY;
    
    alglib::ae_int_t orgSize=newSizeTmp;
    alglib::ae_int_t newSize=orgSizeTmp;
    inputX.setlength(orgSize);
    inputY.setlength(orgSize);
    outputX.setlength(newSize);
    outputY.setlength(newSize);
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        inputX[i]=outputXTmp[i];
        inputY[i]=outputYTmp[i]+delta[i];
    }
    
    for(label i=0; i<NCornerPoints; i++)
    {
        outputX[i]=inputXTmp[i];
    }
    
    alglib::spline1dconvcubic
    (
        inputX,
        inputY,
        orgSize,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputX,
        newSize,
        outputY
    );
    
    for(label i=0; i<NCornerPoints; i++)
    {
        locBcPointsNew[i].component(vector::X)=outputX[i];
        locBcPointsNew[i].component(vector::Y)=outputY[i];
        locBcPointsNew[i].component(vector::Z)=0;
    }
    
    //-Transform from local back to global, take care of the index mismatch.
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            bcPointsNew[i]=twoDLocToGloTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                locBcPointsNew[marker[i]]
            );
            //bcPointsNew[i].component(vector::Z)=0;
            
            /*
            bcPointsNew[markerA[i]].component(vector::X)=-bcPointsNew[i].component(vector::X);
            bcPointsNew[markerA[i]].component(vector::Y)= bcPointsNew[i].component(vector::Y);
            bcPointsNew[markerB[i]].component(vector::X)=-bcPointsNew[i].component(vector::X);
            bcPointsNew[markerB[i]].component(vector::Y)=-bcPointsNew[i].component(vector::Y);
            bcPointsNew[markerC[i]].component(vector::X)= bcPointsNew[i].component(vector::X);
            bcPointsNew[markerC[i]].component(vector::Y)=-bcPointsNew[i].component(vector::Y);
            bcPointsNew[markerD[i]].component(vector::X)= bcPointsNew[i].component(vector::X);
            bcPointsNew[markerD[i]].component(vector::Y)= bcPointsNew[i].component(vector::Y);
            bcPointsNew[markerE[i]].component(vector::X)=-bcPointsNew[i].component(vector::X);
            bcPointsNew[markerE[i]].component(vector::Y)= bcPointsNew[i].component(vector::Y);
            bcPointsNew[markerF[i]].component(vector::X)=-bcPointsNew[i].component(vector::X);
            bcPointsNew[markerF[i]].component(vector::Y)=-bcPointsNew[i].component(vector::Y);
            bcPointsNew[markerG[i]].component(vector::X)= bcPointsNew[i].component(vector::X);
            bcPointsNew[markerG[i]].component(vector::Y)=-bcPointsNew[i].component(vector::Y);
            */
            
            //Info<< "bcPointsNew[" << i << "]= " << tab << bcPointsNew[i] << tab << "|" << tab
            //    << "locBcPointsNew[" << marker[i] << "]= " << tab << locBcPointsNew[marker[i]] << endl;
            /*
            printf
            (
                "New[%d]=(%+2.5f,%+2.5f,%+2.5f) | "
                "Org[%d]=(%+2.5f,%+2.5f,%+2.5f) | "
                "LocNew[%2d]=(%+2.5f,%+2.5f,%+2.5f) | "
                "LocOrg[%2d]=(%+2.5f,%+2.5f,%+2.5f)\n",
                i,
                bcPointsNew[i].component(vector::X),
                bcPointsNew[i].component(vector::Y),
                bcPointsNew[i].component(vector::Z),
                i,
                bcPointsOrg[i].component(vector::X),
                bcPointsOrg[i].component(vector::Y),
                bcPointsOrg[i].component(vector::Z),
                marker[i],
                locBcPointsNew[marker[i]].component(vector::X),
                locBcPointsNew[marker[i]].component(vector::Y),
                locBcPointsNew[marker[i]].component(vector::Z),
                marker[i],
                locBcPointsOrg[marker[i]].component(vector::X),
                locBcPointsOrg[marker[i]].component(vector::Y),
                locBcPointsOrg[marker[i]].component(vector::Z)
            );
            */
        }
    }
    //Info<<endl;
    
//------------------------------------------------------------------------------
// Step (5) Working on the 2nd corner
//------------------------------------------------------------------------------

    mm=2;
    nn=3;
    
    ctrlPointsTmp[0].component(vector::X)=
        -0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[0].component(vector::Y)=0;
    ctrlPointsTmp[0].component(vector::Z)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::X)=
        0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Y)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Z)=0;
    for (label i=1; i<=ctrlPoints_.size(); i++)
    {
        ctrlPointsTmp[i].component(vector::X)=
            ctrlPoints_[i-1].component(vector::X);
        ctrlPointsTmp[i].component(vector::Y)=
            ctrlPoints_[i-1].component(vector::Y);
        ctrlPointsTmp[i].component(vector::Z)=
            ctrlPoints_[i-1].component(vector::Z);
    }
    pointField ctrlPointsTmpp=ctrlPointsTmp;
    for  (label i=0; i<ctrlPointsTmp.size(); i++)
    {
        ctrlPointsTmp[i].component(vector::X)=
            ctrlPointsTmpp[ctrlPointsTmp.size()-i-1].component(vector::X);
        ctrlPointsTmp[i].component(vector::Y)=
            ctrlPointsTmpp[ctrlPointsTmp.size()-i-1].component(vector::Y);
        ctrlPointsTmp[i].component(vector::Z)=
            ctrlPointsTmpp[ctrlPointsTmp.size()-i-1].component(vector::Z);
    }
    
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) > lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) < lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            NCornerPoints=NCornerPoints+1;
        }
    }
    
    marker=0;
    markerInverse=0;
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
            
        )
        {
            marker[i] = NCornerPoints;
            markerInverse[NCornerPoints]=i;
            
            locBcPointsOrg[NCornerPoints]=twoDGloToLocTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                bcPointsNew[i]
            );
            NCornerPoints=NCornerPoints+1;
        }
    }
    
    
    //-Call ALGLIB to do the job!
    
    for(label i=0; i<NCornerPoints; i++)
    {
        inputXTmp[i]=locBcPointsOrg[i].component(vector::X);
        inputYTmp[i]=locBcPointsOrg[i].component(vector::Y);
    }
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        outputXTmp[i]=ctrlPointsTmp[i].component(vector::X);
    }
        
    alglib::spline1dconvcubic
    (
        inputXTmp,
        inputYTmp,
        orgSizeTmp,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputXTmp,
        newSizeTmp,
        outputYTmp
    );
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        inputX[i]=outputXTmp[i];
        inputY[i]=outputYTmp[i]+delta[i];
    }
    
    for(label i=0; i<NCornerPoints; i++)
    {
        outputX[i]=inputXTmp[i];
    }
    
    alglib::spline1dconvcubic
    (
        inputX,
        inputY,
        orgSize,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputX,
        newSize,
        outputY
    );
    
    for(label i=0; i<NCornerPoints; i++)
    {
        locBcPointsNew[i].component(vector::X)=outputX[i];
        locBcPointsNew[i].component(vector::Y)=outputY[i];
    }
    
    //-Transform from local back to global, take care of the index mismatch.
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            bcPointsNew[i]=twoDLocToGloTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                locBcPointsNew[marker[i]]
            );
        }
    }
            Info<< "lockEndAngles[mm]" << lockEndAngles[mm] << endl;
            Info<< "lockEndAngles[nn]" << lockEndAngles[nn] << endl;

//------------------------------------------------------------------------------
// Step (6) Working on the 3rd corner
//------------------------------------------------------------------------------

    mm=4;
    nn=5;
    
    ctrlPointsTmp[0].component(vector::X)=
        -0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[0].component(vector::Y)=0;
    ctrlPointsTmp[0].component(vector::Z)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::X)=
        0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Y)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Z)=0;
    for (label i=1; i<=ctrlPoints_.size(); i++)
    {
        ctrlPointsTmp[i].component(vector::X)=
            ctrlPoints_[i-1].component(vector::X);
        ctrlPointsTmp[i].component(vector::Y)=
            ctrlPoints_[i-1].component(vector::Y);
        ctrlPointsTmp[i].component(vector::Z)=
            ctrlPoints_[i-1].component(vector::Z);
    }
    
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) > lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) < lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            NCornerPoints=NCornerPoints+1;
        }
    }
    
    marker=0;
    markerInverse=0;
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            marker[i] = NCornerPoints;
            markerInverse[NCornerPoints]=i;
            
            locBcPointsOrg[NCornerPoints]=twoDGloToLocTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                bcPointsNew[i]
            );
            NCornerPoints=NCornerPoints+1;
        }
    }
    
    
    //-Call ALGLIB to do the job!
    
    for(label i=0; i<NCornerPoints; i++)
    {
        inputXTmp[i]=locBcPointsOrg[i].component(vector::X);
        inputYTmp[i]=locBcPointsOrg[i].component(vector::Y);
    }
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        outputXTmp[i]=ctrlPointsTmp[i].component(vector::X);
    }
        
    alglib::spline1dconvcubic
    (
        inputXTmp,
        inputYTmp,
        orgSizeTmp,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputXTmp,
        newSizeTmp,
        outputYTmp
    );
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        inputX[i]=outputXTmp[i];
        inputY[i]=outputYTmp[i]+delta[i];
    }
    
    for(label i=0; i<NCornerPoints; i++)
    {
        outputX[i]=inputXTmp[i];
    }
    
    alglib::spline1dconvcubic
    (
        inputX,
        inputY,
        orgSize,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputX,
        newSize,
        outputY
    );
    
    for(label i=0; i<NCornerPoints; i++)
    {
        locBcPointsNew[i].component(vector::X)=outputX[i];
        locBcPointsNew[i].component(vector::Y)=outputY[i];
    }
    
    //-Transform from local back to global, take care of the index mismatch.
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            bcPointsNew[i]=twoDLocToGloTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                locBcPointsNew[marker[i]]
            );
        }
    }

//------------------------------------------------------------------------------
// Step (7) Working on the 4th corner
//------------------------------------------------------------------------------

    mm=6;
    nn=7;
    
    ctrlPointsTmp[0].component(vector::X)=
        -0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[0].component(vector::Y)=0;
    ctrlPointsTmp[0].component(vector::Z)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::X)=
        0.5*mag(lockEndPoints_[mm]-lockEndPoints_[nn]);
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Y)=0;
    ctrlPointsTmp[ctrlPoints_.size()+1].component(vector::Z)=0;
    for (label i=1; i<=ctrlPoints_.size(); i++)
    {
        ctrlPointsTmp[i].component(vector::X)=
            ctrlPoints_[i-1].component(vector::X);
        ctrlPointsTmp[i].component(vector::Y)=
            ctrlPoints_[i-1].component(vector::Y);
        ctrlPointsTmp[i].component(vector::Z)=
            ctrlPoints_[i-1].component(vector::Z);
    }
    ctrlPointsTmpp=ctrlPointsTmp;
    for  (label i=0; i<ctrlPointsTmp.size(); i++)
    {
        ctrlPointsTmp[i].component(vector::X)=
            ctrlPointsTmpp[ctrlPointsTmp.size()-i-1].component(vector::X);
        ctrlPointsTmp[i].component(vector::Y)=
            ctrlPointsTmpp[ctrlPointsTmp.size()-i-1].component(vector::Y);
        ctrlPointsTmp[i].component(vector::Z)=
            ctrlPointsTmpp[ctrlPointsTmp.size()-i-1].component(vector::Z);
    }
    
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) > lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) < lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            NCornerPoints=NCornerPoints+1;
        }
    }
    
    marker=0;
    markerInverse=0;
    NCornerPoints=0;
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            marker[i] = NCornerPoints;
            markerInverse[NCornerPoints]=i;
            
            locBcPointsOrg[NCornerPoints]=twoDGloToLocTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                bcPointsNew[i]
            );
            NCornerPoints=NCornerPoints+1;
        }
    }
    
    
    //-Call ALGLIB to do the job!
    
    for(label i=0; i<NCornerPoints; i++)
    {
        inputXTmp[i]=locBcPointsOrg[i].component(vector::X);
        inputYTmp[i]=locBcPointsOrg[i].component(vector::Y);
    }
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        outputXTmp[i]=ctrlPointsTmp[i].component(vector::X);
    }
        
    alglib::spline1dconvcubic
    (
        inputXTmp,
        inputYTmp,
        orgSizeTmp,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputXTmp,
        newSizeTmp,
        outputYTmp
    );
    
    for(label i=0; i<ctrlPointsTmp.size(); i++)
    {
        inputX[i]=outputXTmp[i];
        inputY[i]=outputYTmp[i]+delta[i];
    }
    
    for(label i=0; i<NCornerPoints; i++)
    {
        outputX[i]=inputXTmp[i];
    }
    
    alglib::spline1dconvcubic
    (
        inputX,
        inputY,
        orgSize,
        boundType,
        boundTypeLeftValue,
        boundType,
        boundTypeRightValue,
        outputX,
        newSize,
        outputY
    );
    
    for(label i=0; i<NCornerPoints; i++)
    {
        locBcPointsNew[i].component(vector::X)=outputX[i];
        locBcPointsNew[i].component(vector::Y)=outputY[i];
    }
    
    //-Transform from local back to global, take care of the index mismatch.
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) >= lockEndAngles[mm] &&
            pointGlobalAngle(globalOrigin_,bcPointsNew[i]) <= lockEndAngles[nn] &&
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            bcPointsNew[i]=twoDLocToGloTransSpecial
            (
                globalOrigin_,
                globalYUnit_,
                lockEndPoints_[nn],
                lockEndPoints_[mm],
                locBcPointsNew[marker[i]]
            );
        }
    }
    
    
    // 3D considerations
    forAll (bcPtsIDsCur, i)
    {
        if
        (
            bcPointsNew[i].component(vector::Z) == 0
        )
        {
            forAll (bcPtsIDsCur, j)
            {
                if
                (
                    bcPointsOrg[j].component(vector::X)==
                        bcPointsOrg[i].component(vector::X) &&
                    bcPointsOrg[j].component(vector::Y)==
                        bcPointsOrg[i].component(vector::Y) &&
                    bcPointsOrg[j].component(vector::Z)!=
                        bcPointsOrg[i].component(vector::Z)
                )
                {
                    bcPointsNew[j].component(vector::X)=
                        bcPointsNew[i].component(vector::X);
                    bcPointsNew[j].component(vector::Y)=
                        bcPointsNew[i].component(vector::Y);
                    bcPointsNew[j].component(vector::Z)=
                        bcPointsOrg[i].component(vector::Z);
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
// Step (8) Calculate the absolute displacement for the whole bc
//------------------------------------------------------------------------------
    //Info << "Hello from here over " << endl;
    
    vectorField dispCur(bcPointsOrg.size(),vector::zero);

if
(
    t.value() >= deformStaTime &&
    t.value() < deformEndTime
)
{
    scalar aaa=
        (deformEndTime-deformStaTime-t.deltaTValue())/t.deltaTValue();
    scalar bbb=
        (t.value()-deformStaTime)/t.deltaTValue();
    dispCur = (bcPointsNew - bcPointsOrg)*bbb/aaa;
    Info<< "Finished " << bbb << " steps in " << aaa << " steps" << endl;
}
else
{
    dispCur = bcPointsCur - bcPointsOrg;
}
    
    forAll (bcPtsIDsOrg, i)
    {
        dispCur[i].component(vector::Z)=0;
    }
    
    Field<vector>::operator=(dispCur);

    fixedValuePointPatchVectorField::updateCoeffs();
}


void udfTypeCMotionWallPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    initialPoints_.writeEntry("initialPoints", os);
    lockEndPoints_.writeEntry("lockEndPoints", os);
    ctrlPoints_.writeEntry("ctrlPoints", os);
    os.writeKeyword("staChangingTime")
        << staChangingTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("endChangingTime")
        << endChangingTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("globalOrigin")
        << globalOrigin_ << token::END_STATEMENT << nl;
    os.writeKeyword("globalYUnit")
        << globalYUnit_ << token::END_STATEMENT << nl;
    os.writeKeyword("endBCType")
        << endBCType_ << token::END_STATEMENT << nl;
    os.writeKeyword("endBCTypeValL")
        << endBCTypeValL_ << token::END_STATEMENT << nl;
    os.writeKeyword("endBCTypeValR")
        << endBCTypeValR_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    udfTypeCMotionWallPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

