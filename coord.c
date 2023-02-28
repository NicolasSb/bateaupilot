#include "coord.h"


double degToRad(double deg)
{
	return deg * M_PI / 180;
}

double radToDeg(double rad)
{
	return rad * 180 / M_PI;
}

double geocentricLatitude(double lat)
{
    // Convert geodetic latitude 'lat' to a geocentric latitude 'clat'.
    // Geodetic latitude is the latitude as given by GPS.
    // Geocentric latitude is the angle measured from center of Earth between a point and the equator.
    double clat = atan((1.0 - ECCENTRICITY_SQUARED) * tan(lat));
    return clat;
}


double earthRadiusInMeters(double latitudeRadians)
{
    // latitudeRadians is geodetic, i.e. that reported by GPS.
    double c = cos(latitudeRadians);
    double s = sin(latitudeRadians);
    double t1 = c * RADIUS * RADIUS;
    double t2 = RADIUS_POLE * RADIUS_POLE * s;
    double t3 = RADIUS * c;
    double t4 = RADIUS_POLE * s;
    return sqrt((t1 * t1 + t2 * t2) / (t3 * t3 + t4 * t4));
}


void locationToPoint(Coord* c, Point* xyzpoint, Vector* norm, double * radius)
{
    // Convert (lat, lon, elv) to (x, y, z).
    double lat = degToRad(c->latitude);
    double lon = degToRad(c->longitude);
    *radius = earthRadiusInMeters(lat);
    double clat = geocentricLatitude(lat);

    double cosLon = cos(lon);
    double sinLon = sin(lon);
    double cosLat = cos(clat);
    double sinLat = sin(clat);
    xyzpoint->x = *radius * cosLon * cosLat;
    xyzpoint->y = *radius * sinLon * cosLat;
    xyzpoint->z = *radius * sinLat;

    // We used geocentric latitude to calculate (x,y,z) on the Earth's ellipsoid.
    // Now we use geodetic latitude to calculate normal vector from the surface, to correct for elevation.
    double cosGlat = cos(lat);
    double sinGlat = sin(lat);

    norm->x = cosGlat * cosLon;
    norm->y = cosGlat * sinLon;
    norm->z = sinGlat;

    xyzpoint->x += c->elevation * norm->x;
    xyzpoint->y += c->elevation * norm->y;
    xyzpoint->z += c->elevation * norm->z;
}



void pointToLocation(Coord* c, Point* xyzpoint, double* radius)
{
    // Convert (x, y, z) to (lat, lon, elv).
    double x = xyzpoint->x;
    double y = xyzpoint->y;
    double z = xyzpoint->z;
    double h = c->elevation;
    double p = sqrt(x * x + y * y);
    double lat = atan((z/p) / (1 - ECCENTRICITY_SQUARED * (*radius / (*radius + h))));
    double lon = atan2(y, x); 
    c->latitude = radToDeg(lat);
    c->longitude = radToDeg(lon);
    c->elevation = h;
}

void rotateGlobe(Coord* b, Coord* a, double bradius, Point* brotated, double* brad)
{
    // Get modified coordinates of 'b' by rotating the globe so that 'a' is at lat=0, lon=0.
    Coord br = { b->latitude, (b->longitude - a->longitude), b->elevation };
    Point xyzpoint;
    Vector norm;
    double radius;
    locationToPoint(&br, &xyzpoint, &norm, &radius);

    // Rotate brp cartesian coordinates around the z-axis by a.lon degrees,
    // then around the y-axis by a.lat degrees.
    // Though we are decreasing by a.lat degrees, as seen above the y-axis,
    // this is a positive (counterclockwise) rotation (if B's longitude is east of A's).
    // However, from this point of view the x-axis is pointing left.
    // So we will look the other way making the x-axis pointing right, the z-axis
    // pointing up, and the rotation treated as negative.

    double alat = geocentricLatitude(degToRad(-a->latitude));
    double acos = cos(alat);
    double asin = sin(alat);

    brotated->x = (xyzpoint.x * acos) - (xyzpoint.z * asin);
    brotated->y = xyzpoint.y;
    brotated->z = (xyzpoint.x * asin) + (xyzpoint.z * acos);

    *brad = bradius;
}

void undoRotation(Coord* a, Point* brotated, Point * unrotated)
{
    double brx = brotated->x;
    double bry = brotated->y;
    double brz = brotated->z;

    double alat = geocentricLatitude(degToRad(-a->latitude)); //maybe minus goes away
    double acos = cos(alat);
    double asin = sin(alat);

    unrotated->y = bry;
    double den = acos + ((asin * asin) / (acos));
    double num = brz - (asin / acos) * brx;
    unrotated->z = num / den;
    unrotated->x = (brx + asin * unrotated->z) / acos;
}


int normalizeVectorDiff(Point b, Point a, Vector* normalized, double* radius)
{
    // Calculate norm(b-a), where norm divides a vector by its length to produce a unit vector.
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double dz = b.z - a.z;
    double dist2 = dx * dx + dy * dy + dz * dz;
    if (dist2 == 0) {
        return 1;
    }
    double dist = sqrt(dist2);
    normalized->x = dx / dist;
    normalized->y = dy / dist;
    normalized->z = dz / dist;

    *radius = 1.0;

    return 0;
}

double distance(Point ap, Point bp)
{
    double dx = ap.x - bp.x;
    double dy = ap.y - bp.y;
    double dz = ap.z - bp.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

void calculate(Coord a, Coord b, double * dist, double * azimuth, double * altitude)
{
    Point coordA, coordB, brotated;
    Vector normA, normB, normalizedDiff;
    double radA, radB, radBrotated, radDiff = 0;
    locationToPoint(&a, &coordA, &normA, &radA);
    locationToPoint(&b, &coordB, &normB, &radB);
    *dist = 0.001 * distance(coordA, coordB);

    // Let's use a trick to calculate azimuth:
    // Rotate the globe so that point A looks like latitude 0, longitude 0.
    // We keep the actual radii calculated based on the oblate geoid,
    // but use angles based on subtraction.
    // Point A will be at x=radius, y=0, z=0.
    // Vector difference B-A will have dz = N/S component, dy = E/W component.
    rotateGlobe(&b, &a, radB, &brotated, &radBrotated);
    if (brotated.z * brotated.z + brotated.y * brotated.y > 1.0e-6)
    {
        double theta = radToDeg(atan2(brotated.z, brotated.y));
        *azimuth = 90.0 - theta;
        if (*azimuth < 0.0) {
            *azimuth += 360.0;
        }
        if (*azimuth > 360.0) {
            *azimuth -= 360.0;
        }

        int err = normalizeVectorDiff(coordA, coordB, &normalizedDiff, &radDiff);
        if (err == 0)
        {
            // Calculate altitude, which is the angle above the horizon of B as seen from A.
            // Almost always, B will actually be below the horizon, so the altitude will be negative.
            // The dot product of bma and norm = cos(zenith_angle), and zenith_angle = (90 deg) - altitude.
            // So altitude = 90 - acos(dotprod).
            *altitude = 90.0 - (180.0 / M_PI) * acos(normalizedDiff.x * normA.x +
                normalizedDiff.y * normA.y + normalizedDiff.z * normA.z);
        }
    }
}

void goForward(Coord position, Coord* newposition, double azimuth, double speed, double time)
{
    double distance = speed / 3.6 * time;

    double gamma = tan(degToRad(90 - azimuth));

    double radius;
    Point sourcePos;
    Vector norm;
    locationToPoint(&position, &sourcePos, &norm, &radius);
    int mult = azimuth > 180 ? -1 : 1;
    double by = mult * distance / sqrt(1 + (gamma * gamma));
    double alat = geocentricLatitude(degToRad(-position.latitude));
    double acos = cos(alat);
    double asin = sin(alat);
    double bz = (gamma * by - sourcePos.x*asin) / (acos);
    Point location = { sourcePos.x, by, bz };

    newposition->elevation = position.elevation;
    pointToLocation(newposition, &location, &radius);
    //on réinjecte la longitude d'origine, supprimée pour la rotation
    newposition->longitude += position.longitude;

    //printf("Coordonnees du point avance %lf, %lf\n", newposition->latitude, 
    //     newposition->longitude);
}