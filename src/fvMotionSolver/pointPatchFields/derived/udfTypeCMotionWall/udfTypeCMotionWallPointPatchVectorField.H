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

Class
    Foam::udfTypeCMotionWallPointPatchVectorField

Description
    Given a rectangular cylinder (2D), width is B, height is H. B is not 
    necessary to be larger than H.
    So the initial shape is a rectangular cylinder (2D), with its round corner 
    radius (RCR) is zero. Now we increase this RCR, until it reaches to its 
    maxima, that is (H/2), suppose B>H. 
    To keep the maximum conformity, angles are first calculated from the old
    frame, and then projected this point with the same angle to the new curve.
    This might be not the best way, but surely this is the easiest way so far.
    And quiete staightforward in our case. DW, 9/Nov/2011
    
    Modification (DW, Nov/28/2011)
    The initial shape is not necessarily to be rectangular, but it must be
    circumscribed by a rectangular, so to make the c++ realization simpler.
    
    Modification (DW, Dec/01/2011)
    Now it can handle the four corners, so it works finally.
    The initial mesh should be carefully built, so to make more sense, 
    especially around the corners.
    
TO-DO
    Given a arbitrary shape and arbitrary mesh.
    
Some Remarks backup
    //-These are the lengths from the left/right lock ends to the rec corner.
    //-Note, the idea is originally based on a rectangular shape with sharp 
    // corner, however the initial shape is not necessary to be a sharp corner,
    // any kind of shape is okay, as long as their two end tangential edges
    // can be extended to converge on a perpendicular sharp corner.
    //scalar lockEndLL = mag(ctrlPoints_[NS-1].component(vector::X)) * scaleUp;
    //scalar lockEndLR = mag(ctrlPoints_[NE-1].component(vector::X)) * scaleUp;
    
    //-The prerequisite here is important, that we need the control points be 
    // only able to move in the local y derection.
    //-Note the local coordinate is defined as to devide the shape corner to
    // two 45-degree angles, and then the deviding line as the y axis, and 
    // perpendicular to it is the x axis, so to make the local x=0 line goes
    // through the "sharp" corner.
    //-So, ctrlPtsLocalX and delta are required, but ctrlPtsLocalY is not!
    //-The user are requested to input the local x position, and its 
    // corresponding delta values.
    //-And note, (1) Local x is NOT the distance between the points. (2) local 
    // x=0 is the line goes through the sharp corner.
    //
    //-So far, I would wish that these control points be evenly distributed and
    // starting and ending with the lock ends. The lock ends, means its delta 
    // is zero. DW, Nov/28/2011

SourceFiles
    udfTypeCMotionWallPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef udfTypeCMotionWallPointPatchVectorField_H
#define udfTypeCMotionWallPointPatchVectorField_H

#include "fixedValuePointPatchFields.H"
#include "interpolation.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
  Class udfTypeCMotionWallPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class udfTypeCMotionWallPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        //- Initial positions of points on the patch
        pointField initialPoints_;

        //- Currently, only works for quadrilateral shape, with four corners.
        //  Please input the exact location of eight lock-ends, its order
        //  should be input cloclwisely, i.e. SEL-SER-SWL-SWR-NWL-NWR-NEL-NER.
        //  TO-DO, polygon handling.
        pointField lockEndPoints_;

        //- The control points or ends of each segments
        //  The input of ctrlPoints much be in local coordinate, see manual.
        //  The input order must from left to right.
        //  Should not include the two lock end points.
        //  The y() component in this ctrlPoints_ is actually not the real y,
        //  but the target delta displacement.
        pointField ctrlPoints_;

        //- 
        scalar staChangingTime_;

        //- 
        scalar endChangingTime_;
        
        //- 
        point globalOrigin_;
        
        //- 
        vector globalYUnit_;

        //- This is used in ALGLIB.
        scalar endBCType_;

        //- This is used in ALGLIB.
        scalar endBCTypeValL_;

        //- This is used in ALGLIB.
        scalar endBCTypeValR_;


    // Private Member Functions

        //- I do not use const point&, since there is a warning there.
        const point twoDGloToLocTransSpecial
        (
            const point& globalOrigin,
            const vector& globalYUnit,
            const point& lockEndLt,
            const point& lockEndRt,
            const point& globalPoint
        ) const;

        //- 
        const point twoDLocToGloTransSpecial
        (
            const point& globalOrigin,
            const vector& globalYUnit,
            const point& lockEndLt,
            const point& lockEndRt,
            const point& globalPoint
        ) const;

        //- 
        // I removed the const ahead of this here, so to remove the warning.
        // I also deleted the ap.h and ap.cpp files, line-709 and line-710 const 
        // to remove the warning.
        scalar pointGlobalAngle
        (
            const point& globalOrigin,
            const point& globalPoint
        ) const;

public:

    //- Runtime type information
    TypeName("udfTypeCMotionWall");


    // Constructors

        //- Construct from patch and internal field
        udfTypeCMotionWallPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        udfTypeCMotionWallPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        udfTypeCMotionWallPointPatchVectorField
        (
            const udfTypeCMotionWallPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new udfTypeCMotionWallPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        udfTypeCMotionWallPointPatchVectorField
        (
            const udfTypeCMotionWallPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new udfTypeCMotionWallPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

