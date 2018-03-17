"""
Author: Brian Leeson
this file takes all files in the source directory, parses the lat lon data into a a format that can be transformed into a vtk file.
this file will likely write into something that looks like a C array

possible lat lon corner boundries

SW Corner
43.974561, -123.226593

NE Corner
44.109570, -122.923095

current source folder contains runs for JAN, FEB

"""

import os


def convert(srcPath, destPath, latRange, longRange, eleRange):
    """
    function takes a path to the location of one or more .gpx files.
    treats lat, long, elevation in each file as (x,y,z) coordinates
    writes coords to a C array of the form:
    float *coords = {x1, y1, z1. x2, y2, z2, ...}
    :param srcPath: string path to gpx files
    :param destPath: string path to place C file
    :param latRange: tuple. (float low, float high)
    :param longRange: tuple. (float low, float high)
    :param eleRange: tuple. (float low, float high)
    :return: None
    """

    allFilePoints = []
    # get each gpx file
    for gpxFile in os.listdir(srcPath):
        if(gpxFile.split(".")[1] == "gpx"):
            # extract all points[(x1, y1, z1), ...] out of the file
            filePoints = extractXYZ(gpxFile, srcPath)
            allFilePoints.extend(filePoints)

    # write points within lat, long, ele range to file
    with open(destPath, "w") as filePointer:
        prefix = "static float allPoints[{}] = {}\n".format(len(allFilePoints), "{")
        filePointer.write(prefix)
        for pointIndex in range(len(allFilePoints)):
            lat, long, ele = allFilePoints[pointIndex]
            inLatRange = (latRange[0] < lat) and (lat < latRange[1])
            inLonRange = (longRange[0] < long) and (long < longRange[1])
            inEleRange = (eleRange[0] < ele) and (ele < eleRange[1])

            # if in range, write to file
            if (inLatRange and inLonRange and inEleRange):
                filePointer.write("{}, {}, {}, ".format(lat, long, ele))

                #  add newline every 3rd point.
                if ((pointIndex + 1) % 3 == 0):
                    filePointer.write("\n")

        postfix = "}"
        filePointer.write(postfix)

    return None


def extractXYZ(filename, path):
    """
    takes a filename of a gpx file and the path to that file.
    extracts the lat long and elevation from each line.
    returns a list of tuples of floats of the form
    [(lat1, long1, ele1), (lat2, long2, ele2), ...]
    :param filename:
    :param path:
    :return: list of tuples
    """
    relativePath = path + filename
    with open(relativePath, "r") as filePointer:
        # NOTE: gpx puts all data on one line
        startingTag = "<trkpt "
        dataList = filePointer.read().split(startingTag)

        pointList = []
        # iterate over each point in dataList and create our list of tuples
        for point in dataList[1:]:  # skip the first element, the header
            # extract lat, long, ele
            latStart = point.find("\"") + 1
            latEnd = point[latStart:].find("\"") + latStart
            lat = float(point[latStart: latEnd])

            point = point.split("lon=\"")[1]

            longEnd = point.find("\"")
            long = float(point[:longEnd])

            ele = float(point.split("ele>")[1][:-2])

            pointList.append((lat, long, ele))

        return pointList


def main():
    srcPath = "./sourceGPX/"
    destPath = "./pointArray.c"
    latRange = (43.974561, 44.109570)
    longRange = (-123.226593, -122.923095)
    eleRange = (0.0, 2500)  # tallest thing in eugene = 2,058

    convert(srcPath, destPath, latRange, longRange, eleRange)

    print("Done")
    return None

if __name__ == "__main__":
    main()
