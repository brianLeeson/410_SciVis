/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

int main()
{
    int  i, j;
	float iso_val = 3.2;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
	int numSegments[16] = {0,1,1,1,1,1,2,1,1,2,1,1,1,1,1,0};
	int lup[16][4];
    lup[0][0] = lup[0][1] = lup[0][2] = lup[0][3] = -1; //0
    lup[1][0] = 0; lup[1][1] = 3; lup[1][2] = -1; lup[1][3] = -1; //1
    lup[2][0] = 0; lup[2][1] = 1; lup[2][2] = -1; lup[2][3] = -1; //2
    lup[3][0] = 1; lup[3][1] = 3; lup[3][2] = -1; lup[3][3] = -1; //3
    lup[4][0] = 2; lup[4][1] = 3; lup[4][2] = -1; lup[4][3] = -1; //4
    lup[5][0] = 0; lup[5][1] = 2; lup[5][2] = -1; lup[5][3] = -1; //5
    lup[6][0] = 0; lup[6][1] = 1; lup[6][2] = 2; lup[6][3] = 3; //6
    lup[7][0] = 1; lup[7][1] = 2; lup[7][2] = -1; lup[7][3] = -1; //7
    lup[8][0] = 1; lup[8][1] = 2; lup[8][2] = -1; lup[8][3] = -1; //8
    lup[9][0] = 0; lup[9][1] = 3; lup[9][2] = 1; lup[9][3] = 2; //9
    lup[10][0] = 0; lup[10][1] = 2; lup[10][2] = -1; lup[10][3] = -1; //10
    lup[11][0] = 2; lup[11][1] = 3; lup[11][2] = -1; lup[11][3] = -1; //11
    lup[12][0] = 1; lup[12][1] = 3; lup[12][2] = -1; lup[12][3] = -1; //12
    lup[13][0] = 0; lup[13][1] = 1; lup[13][2] = -1; lup[13][3] = -1; //13
    lup[14][0] = 0; lup[14][1] = 3; lup[14][2] = -1; lup[14][3] = -1; //14
    lup[15][0] = lup[15][1] = lup[15][2] = lup[15][3] = -1; //15

	//for each cell
	int cellId = 0;
	int num_cells = GetNumberOfCells(dims);
	int idx[2];

	int ll_logical[2], lr_logical[2], ul_logical[2], ur_logical[2];
	float ll_actual[2], lr_actual[2], ul_actual[2], ur_actual[2];
	float ll_scalar, lr_scalar, ul_scalar, ur_scalar;

	for (int cellId = 0; cellId < num_cells; cellId++)
	{
		//figure out what the logical cell index (x,y of lower left corner) is
		GetLogicalCellIndex(idx, cellId, dims);

		//figure out logical point index (x,y) for all 4 corners
		printf("\ncellID: %d, idx (x,y): (%d, %d)\n", cellId, idx[0], idx[1]);

		ll_logical[0] = idx[0];
		ll_logical[1] = idx[1];

		lr_logical[0] = idx[0] + 1;
		lr_logical[1] = idx[1];

		ul_logical[0] = idx[0];
		ul_logical[1] = idx[1] + 1;

		ur_logical[0] = idx[0] + 1;
		ur_logical[1] = idx[1] + 1;

		//range -10 - 10	
		ll_actual[0] = X[ll_logical[0]];
		ll_actual[1] = Y[ll_logical[1]];

		lr_actual[0] = X[lr_logical[0]];
		lr_actual[1] = Y[lr_logical[1]];

		ul_actual[0] = X[ul_logical[0]];
		ul_actual[1] = Y[ul_logical[1]];

		ur_actual[0] = X[ur_logical[0]];
		ur_actual[1] = Y[ur_logical[1]];

		//get the values in the scalar field for each corner
		ll_scalar = F[GetPointIndex(ll_logical, dims)];
		lr_scalar = F[GetPointIndex(lr_logical, dims)];
		ul_scalar = F[GetPointIndex(ul_logical, dims)];
		ur_scalar = F[GetPointIndex(ur_logical, dims)];

		//figure out what case you are in
		unsigned int icase = 0x00;
		if (ll_scalar > iso_val) {icase |= 0x01;}
		if (lr_scalar > iso_val) {icase |= 0x02;}
		if (ul_scalar > iso_val) {icase |= 0x04;}
		if (ur_scalar > iso_val) {icase |= 0x08;}

		printf("in case: %d \n", icase);
		//printf("sizeof(x): %lu\n", sizeof(0x00));
		
		//construct edges
		int nsegments = numSegments[icase];
		int edge1, edge2, i;
		float pt1[2], pt2[2];
		float idist, sdist;
		float FA, FB, FX, A, B;
		for (i = 0 ; i < nsegments; i++)
		{
			edge1 = lup[icase][2*i];
			// Interpolate position along edge1
			if (edge1 == 0)
			{
				printf("edge1: %d\n", edge1);

				FA = ll_scalar;
				FB = lr_scalar;
				A = ll_actual[0];
				B = lr_actual[0];
				FX = iso_val;
				
				pt1[1] = ll_actual[1];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt1[0] =  A + ((idist/sdist) * (B-A));
			}			
			else if (edge1 == 1)
			{
				printf("edge1: %d\n", edge1);

				FA = lr_scalar;
				FB = ur_scalar;
				A = lr_actual[1];
				B = ur_actual[1];
				FX = iso_val;
				
				pt1[0] = lr_actual[0];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt1[1] =  A + ((idist/sdist) * (B-A));
			}
			else if (edge1 == 2)
			{
				printf("edge1: %d\n", edge1);

				FA = ul_scalar;
				FB = ur_scalar;
				A = ul_actual[0];
				B = ur_actual[0];
				FX = iso_val;

				pt1[1] = ul_actual[1];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt1[0] =  A + ((idist/sdist) * (B-A));
			}
			else if (edge1 == 3)
			{
				printf("edge1: %d\n", edge1);

				FA = ll_scalar;
				FB = ul_scalar;
				A = ll_actual[1];
				B = ul_actual[1];
				FX = iso_val;
				
				pt1[0] = ll_actual[0];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt1[1] =  A + ((idist/sdist) * (B-A));
			}
			else {printf("--- ERROR\n"); exit(1);}

			edge2 = lup[icase][2*i+1];
			// Interpolate position along edge2
						if (edge2 == 0)
			{
				printf("edge2: %d\n", edge2);

				FA = ll_scalar;
				FB = lr_scalar;
				A = ll_actual[0];
				B = lr_actual[0];
				FX = iso_val;
				
				pt2[1] = ll_actual[1];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt2[0] =  A + ((idist/sdist) * (B-A));
			}			
			else if (edge2 == 1)
			{
				printf("edge2: %d\n", edge2);

				FA = lr_scalar;
				FB = ur_scalar;
				A = lr_actual[1];
				B = ur_actual[1];
				FX = iso_val;
				
				pt2[0] = lr_actual[0];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt2[1] =  A + ((idist/sdist) * (B-A));
			}
			else if (edge2 == 2)
			{
				printf("edge2: %d\n", edge2);

				FA = ul_scalar;
				FB = ur_scalar;
				A = ul_actual[0];
				B = ur_actual[0];
				FX = iso_val;

				pt2[1] = ul_actual[1];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt2[0] =  A + ((idist/sdist) * (B-A));
			}
			else if (edge2 == 3)
			{
				printf("edge2: %d\n", edge2);

				FA = ll_scalar;
				FB = ul_scalar;
				A = ll_actual[1];
				B = ul_actual[1];
				FX = iso_val;
				
				pt2[0] = ll_actual[0];

				idist = (FX - FA);
				sdist = (FB - FA);
				pt2[1] =  A + ((idist/sdist) * (B-A));
			}
			else {printf("--- ERROR\n"); exit(1);}

			printf("adding point: p1(%f, %f), p2(%f, %f)\n", pt1[0], pt1[1], pt2[0], pt2[1]);
			sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);

		}
	}
	printf("--- DONE\n");

    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //

    iren->Initialize();
    iren->Start();
	printf("Prog End\n");
    pd->Delete();
	
}
