#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


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


float interpolate(const float scalar_A, const float scalar_B, const float actual_A, const float actual_B, const float actual_X){

	float t = (actual_X - actual_A)/(actual_B - actual_A);
	float result = scalar_A + t * (scalar_B - scalar_A);
	return result;
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
                        const float *X, const float *Y, const float *F)
{
	int x_index, y_index;

	// logical index of points bounding containing cell
	int logical_ll[2] = {-1, -1};
	int logical_lr[2] = {-1, -1};
	int logical_ul[2] = {-1, -1};
	int logical_ur[2] = {-1, -1};

	// values in scalar field of points bounding containing cell
	float scalar_ll, scalar_lr, scalar_ul, scalar_ur, top, bot, dif;

//	find the logical index (x,y) of the four corners of the cell that contains pt.
	for(x_index = 0; x_index < (dims[0]-1); x_index++){
		if ((X[x_index] <= pt[0]) && (pt[0] <= X[x_index + 1])){
			logical_ll[0] = x_index;
			logical_lr[0] = x_index + 1;
			logical_ul[0] = x_index;
			logical_ur[0] = x_index + 1;
			break;
		}
	}
	
	for(y_index = 0; y_index < (dims[1]-1); y_index++){
		if ((Y[y_index] <= pt[1]) && (pt[1] <= Y[y_index + 1])){
			logical_ll[1] = y_index;
			logical_lr[1] = y_index;
			logical_ul[1] = y_index + 1;
			logical_ur[1] = y_index + 1;
			break;
		}
	}
	if(logical_ll[0] == -1 || logical_ll[1] == -1){
		return 0;
	}	

//	find values in scalar field of points bounding containing cell
	
	scalar_ll = F[GetPointIndex(logical_ll, dims)];
	scalar_lr = F[GetPointIndex(logical_lr, dims)];
	scalar_ul = F[GetPointIndex(logical_ul, dims)];
	scalar_ur = F[GetPointIndex(logical_ur, dims)];

// 	interpolate find value of pt inside the cell
	// interp bot
	bot = interpolate(scalar_ll, scalar_lr, X[x_index], X[x_index + 1], pt[0]);

	// interp top
	top = interpolate(scalar_ul, scalar_ur, X[x_index], X[x_index + 1], pt[0]);

	// interp dif
	dif = interpolate(bot, top, Y[y_index], Y[y_index + 1], pt[1]);

    return dif;
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
	//printf("%f\n", F);
	int RGB_MAX = 255;
	int R_MIN = 0;
	int G_MIN = 0;
	int B_MIN = 128;

	float r_val = F * (RGB_MAX - R_MIN) + R_MIN;
	float g_val = F * (RGB_MAX - G_MIN) + G_MIN;
	float b_val = F * (RGB_MAX - B_MIN) + B_MIN;

	RGB[0] = r_val;
	RGB[1] = g_val;
	RGB[2] = b_val;
	//printf("%f\n", r_val);
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F = 0:	(0,0,128) 
//        F = 0.5:	(255,255,255) 
//        F = 1: 	(128, 0, 0)
//        and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
	float r_val = 255.0;
	float g_val = 255.0;
	float b_val = 255.0;

	int R_VALS[] = {0, 255, 128};
	int G_VALS[] = {0, 255, 0};
	int B_VALS[] = {128, 255, 0};	

	float cutoff = .5;

	if(F <= cutoff){
		r_val = (F * 2) * (R_VALS[1] - R_VALS[0]) + R_VALS[0];
		g_val = (F * 2) * (G_VALS[1] - G_VALS[0]) + G_VALS[0];
		b_val = (F * 2) * (B_VALS[1] - B_VALS[0]) + B_VALS[0];
	}

	else if(F > cutoff){
		r_val = ((F - cutoff) / cutoff) * (R_VALS[1] - R_VALS[2]) + R_VALS[2];
		g_val = ((F - cutoff) / cutoff) * (G_VALS[1] - G_VALS[2]) + G_VALS[2];
		b_val = ((F - cutoff) / cutoff) * (B_VALS[1] - B_VALS[2]) + B_VALS[2];
	}

	RGB[0] = r_val;
	RGB[1] = g_val;
	RGB[2] = b_val;
}

//From Slides

/*
hsvToRGB(float hue, float saturation, float value)
{
	if(saturation == 0 ) // achromatic (grey)
	{
		r = g = b = v;
	}
	else
	{
		hue /= 60.f;
		// sector 0 to 5
		i = floor( hue );
		f = hue - i;
		// factorial part of h
		p = value * ( 1 - saturation);
		q = value * ( 1 - saturation * f );
		t = value * ( 1 - saturation * ( 1 - f ) );
		switch( i ):
			case 0: r = v; g = t; b = p;
			case 1: r = q; g = v; b p;
			case 2: r = p; g = v; b = t;
			case 3: r = p; g = q; b = v;
			case 4: r = t; g = p; b = v;
			case 5: r = v; g = p; b = q;
		}
}
*/

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

    for (i = 0 ; i < nx ; i++)
    {
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            pt[0] = -9.0 + 18.0 * (((float)i) / ((float)(nx-1)));
            pt[1] = -9.0 + 18.0 * (((float)j) / ((float)(ny-1)));
            //printf("pt: (%f, %f)\n", pt[0], pt[1]);
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            //printf("val: %f\n", f);
            float MAX_SCAL_VAL = 5.02;
            float MIN_SCAL_VAL = 1.2;
            float normalizedF = (f - MIN_SCAL_VAL)/(MAX_SCAL_VAL - MIN_SCAL_VAL); 
            
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }
    }
    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
