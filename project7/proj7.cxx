/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

Credit: I got the idea for a transform array in buildLogicalPointArray from Laura.

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

#include "TriangleList.h"
#include "tricase.cxx"


//GLOBALS
const float ISO_VAL = 3.2;


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
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
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
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
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
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
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
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
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
    idx[0] = pointId%dims[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
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
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
// helper functions
// ****************************************************************************


//given a logical index for a cell, create array of logical values.
void buildLogicalPointArray(int* logicalPointArray, const int* idx)
{
	int transform[3*8] = {0,0,0, 1,0,0, 0,0,1, 1,0,1, 0,1,0, 1,1,0, 0,1,1, 1,1,1}; // order matters
	int i = 0;
	for (i = 0; i<8; i++)
	{
		logicalPointArray[3*i+0] = idx[0] + transform[3*i+0];
		logicalPointArray[3*i+1] = idx[1] + transform[3*i+1];
		logicalPointArray[3*i+2] = idx[2] + transform[3*i+2];
	}
}


//given a logical index for a cell, create array of scalar values. each element in the array is the scalar value of a vertex. array[0] is vertex 0 and so on.
void buildScalarPointArray(float* scalarPointArray, const int* idx, const float* F, const int* dims)
{
	int logicalIndexArray[3*8]; //3 positions for each of the 8 points
	// buildLogicalPointArray
	buildLogicalPointArray(logicalIndexArray, idx);

	//building scalarArray
	int i;
	int pointIndex = -1;
	int logicalPoint[3];
	for(i = 0; i < 8; i++)
	{
		logicalPoint[0] = logicalIndexArray[3*i+0];
		logicalPoint[1] = logicalIndexArray[3*i+1];
		logicalPoint[2] = logicalIndexArray[3*i+2];
		pointIndex = GetPointIndex(logicalPoint, dims);
		scalarPointArray[i] = F[pointIndex];
	}
	
}

// return the current case that you are in.
int findCase(const int cellID, const float ISO_VAL, const float* F, const int* dims)
{
	// get cells logical cell index
	int idx[3];
	GetLogicalCellIndex(idx, cellID, dims);

	// build array of scalar fields values for each point {sp0, sp1, ..., sp7}
	float scalarPointArray[8];  // scalar values for each of the 8 vertecies
	buildScalarPointArray(scalarPointArray, idx, F, dims);

	//determine case
	unsigned int builder = 0x01;
	unsigned int icase = 0x00;
	int i = 0;
	// for value in scalar_fields_values
	for(i = 0; i<8; i++)
	{
		if(scalarPointArray[i] > ISO_VAL){icase |= builder;}
		builder = builder<<1;
	}
	//printf("found icase: %d\n", icase);
	return icase;

}

void interp(float* resultPoint, const int edge, const int* idx, const float* X, const float* Y, const float* Z, const float* F, const int* dims)
{
	// get vertecies that create edge
	//eTV[2*edge+0] and eTV[2*edge+1] to get verts for that edge
	int edgeToVertecies[12*2] = {0,1, 1,3, 2,3, 0,2, 4,5, 5,7, 6,7, 4,6, 0,4, 1,5, 2,6, 3,7};
	int vertA = edgeToVertecies[2*edge+0];
	int vertB = edgeToVertecies[2*edge+1];

	// get logical xyz for vertA and vertB
	int logicalPointArray[3*8]; //3 positions for each of the 8 points
	buildLogicalPointArray(logicalPointArray, idx);
	
	int logicalPoint1[3], logicalPoint2[3];
	logicalPoint1[0] = logicalPointArray[3*vertA + 0];
	logicalPoint1[1] = logicalPointArray[3*vertA + 1];
	logicalPoint1[2] = logicalPointArray[3*vertA + 2];

	logicalPoint2[0] = logicalPointArray[3*vertB + 0];
	logicalPoint2[1] = logicalPointArray[3*vertB + 1];
	logicalPoint2[2] = logicalPointArray[3*vertB + 2];

	// get actual xyz for those points
	float actualPoint1[3], actualPoint2[3];

	actualPoint1[0] = X[logicalPoint1[0]];
	actualPoint1[1] = Y[logicalPoint1[1]];
	actualPoint1[2] = Z[logicalPoint1[2]];

	actualPoint2[0] = X[logicalPoint2[0]];
	actualPoint2[1] = Y[logicalPoint2[1]];
	actualPoint2[2] = Z[logicalPoint2[2]];

	// interpolate

	// get scalar value for those two points
	float FA = F[GetPointIndex(logicalPoint1, dims)];
	float FB = F[GetPointIndex(logicalPoint2, dims)];

	// figure out what your A and B is
	float edgeToIndex[12] = {0, 2, 0, 2, 0, 2, 0, 2, 1, 1, 1, 1};
	int componentIndex = edgeToIndex[edge];
	float A = actualPoint1[componentIndex];
	float B = actualPoint2[componentIndex];

	// initialize result. only one direction changes, the others remain
	resultPoint[0] = actualPoint1[0];
	resultPoint[1] = actualPoint1[1];
	resultPoint[2] = actualPoint1[2];

	float t = (ISO_VAL - FA)/(FB - FA);

	resultPoint[componentIndex] = A+t*(B-A);
}

int main()
{	
	vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj7.vtk");
    rdr->Update();

	int dims[3];
	vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
	float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

	TriangleList tl;

	// Start doing the project

	// iterate through each each cell
	int numCells = GetNumberOfCells(dims);
	int cellNum = 0;

	// for each cell
	for(cellNum = 0; cellNum < numCells; cellNum++)
	{
		// find out what case you are in
		int icase = findCase(cellNum, ISO_VAL, F, dims);
		
		// find out how many triangles you will need to draw - for each triangle to be drawn:
		int idx[3];
		GetLogicalCellIndex(idx, cellNum, dims);
		int i = 0;
		while(triCase[icase][3*i] != -1)
		{
			int edge1 = triCase[icase][3*i+0];
			int edge2 = triCase[icase][3*i+1];
			int edge3 = triCase[icase][3*i+2];

			// interpolate to find the 3 vertecies of the triangle
			float pt1[3], pt2[3], pt3[3];
			interp(pt1, edge1, idx, X, Y, Z, F, dims);
			interp(pt2, edge2, idx, X, Y, Z, F, dims);
			interp(pt3, edge3, idx, X, Y, Z, F, dims);

			// add each vertex to the triangle list (tl)
			tl.AddTriangle(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);

			i++;
		}
	}

	// Stop doing the project

	vtkPolyData *pd = tl.MakePolyData();

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

    iren->Initialize();
    iren->Start();
    pd->Delete();
	
}
