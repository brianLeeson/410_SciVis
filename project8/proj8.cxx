/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

Credit: ---

=========================================================================*/

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

#include <vtkDataSetWriter.h>

int main()
{	
	
	//read data from file
    vtkDataSetReader *reader = vtkDataSetReader::New();
    reader->SetFileName("proj8.vtk");
	

	//countor filter for render #1
	vtkContourFilter *cf = vtkContourFilter::New();
	cf->SetNumberOfContours(1);
	cf->SetValue(0, 2.5);
	cf->SetValue(1, 5.0);
	cf->SetInputConnection(reader->GetOutputPort());

	// mapper #1
	vtkSmartPointer<vtkPolyDataMapper> inputMapper1 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	inputMapper1->SetInputConnection(cf->GetOutputPort());

	// actor #1
	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	actor1->GetProperty()->SetColor(1.0, 0.8941, 0.7686); // bisque
	actor1->SetMapper(inputMapper1);
 

	// render #2

	// render #3

	// render #4

	//this creates a dummy sphere. delete when using real data.
	vtkSmartPointer<vtkPolyData> STUB_inputPolyData;
	vtkSmartPointer<vtkSphereSource> sphereSource =
      vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetThetaResolution(30);
    sphereSource->SetPhiResolution(15);
    sphereSource->Update();
    STUB_inputPolyData = sphereSource->GetOutput();

	vtkSmartPointer<vtkPolyDataMapper> inputMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();

	inputMapper->SetInputData(STUB_inputPolyData);

	// Create input actor
	vtkSmartPointer<vtkActor> testActor = 
		vtkSmartPointer<vtkActor>::New();
	testActor->GetProperty()->SetColor(1.0, 0.8941, 0.7686); // bisque
	testActor->SetMapper(inputMapper);
 

	// Define viewport ranges
	// (xmin, ymin, xmax, ymax)
	double viewport1[4] = {0.0, 0.0, 0.5, 0.5};
	double viewport2[4] = {0.0, 0.5, 0.5, 1.0};
	double viewport3[4] = {0.5, 0.0, 1.0, 0.5};
	double viewport4[4] = {0.5, 0.5, 1.0, 1.0};

	// make renderers
	vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> ren3 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> ren4 = vtkSmartPointer<vtkRenderer>::New();

	ren1->SetViewport(viewport1);
	ren2->SetViewport(viewport2);
	ren3->SetViewport(viewport3);
	ren4->SetViewport(viewport4);

	ren1->AddActor(actor1); //display the sphere
	ren2->AddActor(testActor); //display the sphere
	ren3->AddActor(testActor); //display the sphere
	ren4->AddActor(testActor); //display the sphere
 
	//Add renderer to renderwindow and render
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(ren1);
	renderWindow->AddRenderer(ren2);
	renderWindow->AddRenderer(ren3);
	renderWindow->AddRenderer(ren4);
	renderWindow->SetSize(800, 800);


	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	
	renderWindow->Render();

	interactor->Start();
	
	return EXIT_SUCCESS;

/*
	// Apply filters
	vtkContourFilter *cf = vtkContourFilter::New();
	cf->SetNumberOfContours(1);
	cf->SetValue(0,3.0);
	cf->SetInputConnection(rdr->GetOutputPort());

	
	ren1->AddActor(testActor);


	// make render window
	vtkSmartPointer<vtkRenderWindow> renWin = 
		vtkSmartPointer<vtkRenderWindow>::New();

	// add renderers to render window
	renWin->AddRenderer(ren1);
	//renWin->AddRenderer(ren2);
	//renWin->AddRenderer(ren3);
	//renWin->AddRenderer(ren4);

	renWin->SetSize(600, 600);

	// make render window interactor - needed?
	vtkSmartPointer<vtkRenderWindowInteractor> iren = 
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);


	// this should show the picture
	iren->Initialize();
	iren->Start();

	return 0;	
	*/
}
