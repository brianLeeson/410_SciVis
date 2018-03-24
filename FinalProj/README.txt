Report:

Work Done:
    * Researched and problem solved VisIt readable file formats
    * Research Strava and Polar APIs
    * Wrote python file format converter
    * Scale transform and plot points to create a 3D image of .gpx data
    * Created a pseudo color plot of the elevation of the route
    * 15-20 hours


Problems:
    * VTK problems. Solved by trying other file formats. Landed on .3D. Simple point data
    * APIs were for creating webapps and poorly suited for manipulating personal data.
    * Attempted to write python script to generate picture. Ran into deep issues with path variables, install locations,
        visit .silo databases and python version conflicts.

Summary:
        I was able to create a pseudocolor plot using the data that I collected and formatted. This involved understanding
    and troubleshooting visit, python, vtk. Most problems where found after reading through visit documentation. Through
    this project I have gained additional knowledge on how visit works, how to visualize data and what the different
     parts that make up the data are.

Further work to be done:
        I hope to expand this project both as something for personal use and something to show to potential employers. I would like
    for there to exist a makefile to run or a webpage to go to such that if gpx data was supplied, a user to visualize
    some simple aspect of their own data.
    This would at least involve:
        resolve python2 path conflicts
            correct visit instillation and inclusion in path
            correct version of python2
            fix visit lib import errors
        Makefile construction
        Virtual Environment
        Further learning python calls to create visit filters
        potentially changing file formats to accommodate addition variables