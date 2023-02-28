#ifndef COORD_H
#define COORD_H

#define _USE_MATH_DEFINES
#include <math.h>

#define RADIUS 6378137.0
#define RADIUS_POLE 6356752.3
#define ECCENTRICITY_SQUARED 0.00669437999014
/**
* \brief Structure to define coordinates
**/
typedef struct Coord
{
	double latitude;
	double longitude;
    double elevation;
}Coord;

typedef struct Point
{
    double x;
    double y;
    double z;
}Point;

typedef struct Vector
{
    double x;
    double y;
    double z;
}Vector;


double degToRad(double);

double radToDeg(double);

double geocentricLatitude(double lat);

double earthRadiusInMeters(double latitudeRadians);

void locationToPoint(Coord* c, Point* xyzpoint, Vector* norm, double * radius);

void pointToLocation(Coord* c, Point* xyzpoint, double* radius);

void rotateGlobe(Coord* b, Coord* a, double bradius, Point* brotated, double* brad);

void undoRotation(Coord* a, Point* brotated, Point* unrotated);

int normalizeVectorDiff(Point b, Point a, Vector* normalized, double* radius);

double distance(Point ap, Point bp);

void calculate(Coord a, Coord b, double * distance, double * azimuth, double * altitude);

void goForward(Coord position, Coord* newposition, double azimuth, double speed, double time);



#endif