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

#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkHedgeHog.h>
#include <vtkMaskPoints.h>
#include <vtkStreamTracer.h>
#include "vtkRungeKutta4.h"
#include "vtkLineSource.h"

int main()
{	
	
	//read data from file
    vtkDataSetReader *reader = vtkDataSetReader::New();
    reader->SetFileName("proj8.vtk");
	reader->Update();

	// get scalar range
	double range[2];
	reader->GetOutput()->GetPointData()->GetScalars()->GetRange(range);
	
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
	inputMapper1->SetScalarRange(range);

	// actor #1 - contourFilter
	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	//actor1->GetProperty()->SetColor(1.0, 0.8941, 0.7686); // bisque
	actor1->SetMapper(inputMapper1);
 

	// render #2 - cutter
	vtkSmartPointer<vtkActor> actor2a = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> actor2b = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> actor2c = vtkSmartPointer<vtkActor>::New();
	
	// Create a plane #1 to cut,
	vtkSmartPointer<vtkPlane> plane1 =
		vtkSmartPointer<vtkPlane>::New();
	plane1->SetOrigin(0,0,0);
	plane1->SetNormal(1,0,0);

	// Create cutter #1
	vtkSmartPointer<vtkCutter> cutter1 =
		vtkSmartPointer<vtkCutter>::New();
	cutter1->SetCutFunction(plane1);
	cutter1->SetInputConnection(reader->GetOutputPort());
	cutter1->Update();

	// Create a plane #2 to cut,
	vtkSmartPointer<vtkPlane> plane2 =
		vtkSmartPointer<vtkPlane>::New();
	plane2->SetOrigin(0,0,0);
	plane2->SetNormal(0,1,0);

	// Create cutter #2
	vtkSmartPointer<vtkCutter> cutter2 =
		vtkSmartPointer<vtkCutter>::New();
	cutter2->SetCutFunction(plane2);
	cutter2->SetInputConnection(reader->GetOutputPort());
	cutter2->Update();

	// Create a plane #3 to cut 
	vtkSmartPointer<vtkPlane> plane3 =
		vtkSmartPointer<vtkPlane>::New();
	plane3->SetOrigin(0,0,0);
	plane3->SetNormal(0,0,1);

	// Create cutter #3
	vtkSmartPointer<vtkCutter> cutter3 =
		vtkSmartPointer<vtkCutter>::New();
	cutter3->SetCutFunction(plane3);
	cutter3->SetInputConnection(reader->GetOutputPort());
	cutter3->Update();

	// add mapper to actor
	vtkSmartPointer<vtkPolyDataMapper> cutterMapper1 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	cutterMapper1->SetInputConnection(cutter1->GetOutputPort());
	// set scalar range for mapper1
	cutterMapper1->SetScalarRange(range);	

	vtkSmartPointer<vtkPolyDataMapper> cutterMapper2 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	cutterMapper2->SetInputConnection(cutter2->GetOutputPort());
	// set scalar range for mapper2
	cutterMapper2->SetScalarRange(range);	

	vtkSmartPointer<vtkPolyDataMapper> cutterMapper3 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	cutterMapper3->SetInputConnection(cutter3->GetOutputPort());
	// set scalar range for mapper2
	cutterMapper3->SetScalarRange(range);	

	actor2a->SetMapper(cutterMapper1);
	actor2b->SetMapper(cutterMapper2);
	actor2c->SetMapper(cutterMapper3);
	
	// render #3

	// access the grad portion of the data
	reader->GetOutput()->GetPointData()->SetActiveAttribute("grad",
		vtkDataSetAttributes::VECTORS);
	reader->Update();

	// add vtkMaskPoints filter. set sample rate
	vtkSmartPointer<vtkMaskPoints> hedgeHogMask =
		vtkSmartPointer<vtkMaskPoints>::New();
	hedgeHogMask->SetOnRatio(11);
	hedgeHogMask->SetInputConnection(reader->GetOutputPort());
	hedgeHogMask->Update();

	vtkSmartPointer<vtkHedgeHog> hedgehog = 
		vtkSmartPointer<vtkHedgeHog>::New();
	hedgehog->SetInputData(hedgeHogMask->GetOutput());
	hedgehog->SetScaleFactor(2);
	hedgehog->Update();

	vtkSmartPointer<vtkPolyDataMapper> hedgeHogMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	hedgeHogMapper->SetInputConnection(hedgehog->GetOutputPort());
	hedgeHogMapper->SetScalarRange(range);

	vtkSmartPointer<vtkActor> hedgeHogActor = vtkSmartPointer<vtkActor>::New();
	hedgeHogActor->SetMapper(hedgeHogMapper);

	// render #4
	vtkLineSource *rake = vtkLineSource::New();
	rake->SetPoint1(-9, 0, 0);
	rake->SetPoint2(9, 0, 0);
	rake->SetResolution(19);
	vtkPolyDataMapper *rakeMapper = vtkPolyDataMapper::New();
	rakeMapper->SetInputConnection(rake->GetOutputPort());
	vtkActor *rakeActor = vtkActor::New();
	rakeActor->SetMapper(rakeMapper);

	vtkRungeKutta4 *integ = vtkRungeKutta4::New();
	vtkStreamTracer *streamTracer = vtkStreamTracer::New();
	streamTracer->SetInputConnection(reader->GetOutputPort());
	streamTracer->SetSourceConnection(rake->GetOutputPort());
	streamTracer->SetIntegrator(integ);
	streamTracer->SetMaximumPropagation(10);

	streamTracer->SetInitialIntegrationStep(0.1);

	streamTracer->SetIntegrationDirectionToForward();

	vtkSmartPointer<vtkPolyDataMapper> streamTracerMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	streamTracerMapper->SetInputConnection(streamTracer->GetOutputPort());
	streamTracerMapper->SetScalarRange(range);
	
	vtkActor *streamActor = vtkActor::New();
	streamActor->SetMapper(streamTracerMapper);

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

	ren1->AddActor(actor1); 
	ren2->AddActor(actor2a); 
	ren2->AddActor(actor2b); 
	ren2->AddActor(actor2c);
	ren3->AddActor(hedgeHogActor); 
	ren4->AddActor(streamActor);
 
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
}
