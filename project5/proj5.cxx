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
	lup[0][0] = lup[0][1] = lup[0][2] = lup[0][3] = -1;
	lup[1][0] = 3; lup[1][1] = 0; lup[1][2] = lup[1][3] = -1;
	lup[2][0] = 0; lup[2][1] = 1; lup[2][2] = lup[2][3] = -1;
	lup[3][0] = 3; lup[3][1] = 1; lup[3][2] = lup[3][3] = -1;
	lup[4][0] = 3; lup[4][1] = 2; lup[4][2] = lup[4][3] = -1;
	lup[5][0] = 1; lup[5][1] = 2; lup[5][2] = lup[5][3] = -1;
	lup[6][0] = 0; lup[6][1] = 1; lup[6][2] = 2; lup[6][3] = 3;
	lup[7][0] = 1; lup[7][1] = 2; lup[7][2] = lup[7][3] = -1;
	lup[8][0] = 1; lup[8][1] = 2; lup[8][2] = lup[8][3] = -1;
	lup[9][0] = 0; lup[9][1] = 1; lup[9][2] = 2; lup[9][3] = 3;
	lup[10][0] = 1; lup[10][1] = 2; lup[10][2] = lup[10][3] = -1;
	lup[11][0] = 3; lup[11][1] = 2; lup[11][2] = lup[11][3] = -1;
	lup[12][0] = 3; lup[12][1] = 1; lup[12][2] = lup[12][3] = -1;
	lup[13][0] = 0; lup[13][1] = 1; lup[13][2] = lup[12][3] = -1;
	lup[14][0] = 3; lup[14][1] = 0; lup[14][2] = lup[14][3] = -1;
	lup[15][0] = lup[15][1] = lup[15][2] = lup[15][3] = -1;

	//for each cell index
	int cellId = 0;
	int num_cells = GetNumberOfCells(dims);
	int idx[2];
	int ll_logical[2] = {-1,-1};
	int lr_logical[2] = {-1,-1};
	int ul_logical[2] = {-1,-1};
	int ur_logical[2] = {-1,-1};

	float ll_scalar, lr_scalar, ul_scalar, ur_scalar;
	int X_MIN = X[0];
	int X_MAX = X[dims[0]];
	float X_WIDTH = (float) X_MAX - X_MIN;
	int Y_MIN = Y[0];
	int Y_MAX = Y[dims[1]];
	float Y_HEIGHT = (float) Y_MAX - Y_MIN;

	for (int cellId = 0; cellId < num_cells; cellId++)
	{
		//figure out what the logical cell index (x,y of lower left corner) is
		GetLogicalCellIndex(idx, cellId, dims);

		//figure out logical point index (x,y) for all 4 corners
		ll_logical[0] = idx[0];
		ll_logical[1] = idx[1];

		lr_logical[0] = idx[0] + 1;
		lr_logical[1] = idx[1];

		ul_logical[0] = idx[0];
		ul_logical[1] = idx[1] + 1;

		ur_logical[0] = idx[0] + 1;
		ur_logical[1] = idx[1] + 1;

		//get the values in the scalar field for each corner
		ll_scalar = F[GetPointIndex(ll_logical, dims)];
		lr_scalar = F[GetPointIndex(lr_logical, dims)];
		ul_scalar = F[GetPointIndex(ul_logical, dims)];
		ur_scalar = F[GetPointIndex(ur_logical, dims)];

		//figure out what case you are in
		unsigned char icase = 0x00;
		if (ll_scalar > iso_val) {icase |= 0x01;}
		if (lr_scalar > iso_val) {icase |= 0x02;}
		if (ul_scalar > iso_val) {icase |= 0x04;}
		if (ur_scalar > iso_val) {icase |= 0x08;}

		//printf("%d \n", icase);
		
		//construct edges
		int nsegments = numSegments[icase];
		int edge1, edge2, i;
		float pt1[2], pt2[2];
		float iso_dist, scal_dist;
		for (i = 0 ; i < nsegments ; i++)
		{
			edge1 = lup[icase][2*i];
			// Interpolate position along edge1
			if (edge1 == 0)
			{
				//printf("edge: %d\n", edge1);
				pt1[1] = (float) ll_logical[1];
				iso_dist = fabs(ll_scalar - iso_val);
				scal_dist = fabs(ll_scalar - lr_scalar);
				pt1[0] = ll_logical[0] + iso_dist/scal_dist;
			}			
			else if (edge1 == 1)
			{
				//printf("edge: %d\n", edge1);
				pt1[0] = (float) lr_logical[0];
				iso_dist = fabs(lr_scalar - iso_val);
				scal_dist = fabs(lr_scalar - ur_scalar);
				pt1[1] = lr_logical[1] + iso_dist/scal_dist;
			}
			else if (edge1 == 2)
			{
				//printf("edge: %d\n", edge1);
				pt1[1] = (float) ul_logical[1];
				iso_dist = fabs(ul_scalar - iso_val);
				scal_dist = fabs(ul_scalar - ur_scalar);
				pt1[0] = ul_logical[0] + iso_dist/scal_dist;
			}
			else if (edge1 == 3)
			{
				//printf("edge: %d\n", edge1);
				pt1[0] = (float) ll_logical[0];
				iso_dist = fabs(ll_scalar - iso_val);
				scal_dist = fabs(ll_scalar - ul_scalar);
				pt1[1] = (float) ll_logical[1] + iso_dist/scal_dist;
				//printf("iso_d: %f, scal_d: %f\n", iso_dist, scal_dist);
				//printf("ll_scalar: %f\n", ll_scalar);
				//printf("iso_val: %f\n", iso_val);
				//printf("ll_s - iso_cal = %f\n", ll_scalar - iso_val);
			}
			else printf("--- ERROR\n");


			edge2 = lup[icase][2*i+1];
			// Interpolate position along edge2
			if ( edge2 == 0)
			{
				pt2[1] = (float) ll_logical[1];
				iso_dist = fabs(ll_scalar - iso_val);
				scal_dist = fabs(ll_scalar - lr_scalar);
				pt2[0] = ll_logical[0] + iso_dist/scal_dist;
			}			
			else if ( edge2 == 1)
			{
				pt2[0] = (float) lr_logical[0];
				iso_dist = fabs(lr_scalar - iso_val);
				scal_dist = fabs(lr_scalar - ur_scalar);
				pt2[1] = lr_logical[1] + iso_dist/scal_dist;
			}
			else if ( edge2 == 2)
			{
				pt2[1] = (float) ul_logical[1];
				iso_dist = fabs(ul_scalar - iso_val);
				scal_dist = fabs(ul_scalar - ur_scalar);
				pt2[0] = ul_logical[0] + iso_dist/scal_dist;
			}
			else if ( edge2 == 3)
			{
				pt2[0] = (float) ll_logical[0];
				iso_dist = fabs(ll_scalar - iso_val);
				scal_dist = fabs(ll_scalar - ul_scalar);
				pt2[1] = ll_logical[1] + iso_dist/scal_dist;
			}
			else printf("--- ERROR\n");

			//normalize
			pt1[0] = -10.0 + (pt1[0]/(X_WIDTH)) * 20.0;
			pt1[1] = -10.0 + (pt1[1]/(Y_HEIGHT)) * 20.0;
			pt2[0] = -10.0 + (pt2[0]/(X_WIDTH)) * 20.0;
			pt2[1] = -10.0 + (pt2[1]/(Y_HEIGHT)) * 20.0;
			printf("X_WIDTH: %f, Y_HEIGHT: %f\n", X_WIDTH, Y_HEIGHT);
			// TODO I think the problem is in how pt is calculated initially.
			// normalization probably is correct
			

			//AddLineSegmentToOutput(pt1, pt2);
			printf("point: p1(%f, %f), p2(%f, %f)\n", pt1[0], pt1[1], pt2[0], pt2[1]);
			sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);
		}

		
	}
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

    pd->Delete();
}
