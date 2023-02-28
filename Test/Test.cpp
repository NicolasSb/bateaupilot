#include "pch.h"
#include "CppUnitTest.h"
#include <math.h>
#include "../coord.h"
#include "../coord.c"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Test
{
	TEST_CLASS(Test)
	{
	public:
		
		TEST_METHOD(testdegToRad)
		{
			Assert::IsTrue(degToRad(180) < 3.1416 && degToRad(180)>3.141592);
			Assert::AreEqual(degToRad(0),0.0);
			Assert::IsTrue(degToRad(30) < 0.5236 && degToRad(30) > 0.52359);
		}

		TEST_METHOD(testradtodeg)
		{
			Assert::AreEqual(radToDeg(M_PI), 180.0);
			Assert::AreEqual(radToDeg(0), 0.0);
			Assert::IsTrue(radToDeg(0.52359) < 31 && radToDeg(0.52359) > 29);
		}

		TEST_METHOD(testgeolat)
		{
			for (int i = -90; i < 90; i+=10)
				Assert::AreEqual(geocentricLatitude(degToRad(i)), degToRad(i), 0.01);
		}

		TEST_METHOD(testearthRadiusInMeters)
		{
			Assert::AreEqual(earthRadiusInMeters(degToRad(0)), RADIUS, 10);
			Assert::AreEqual(earthRadiusInMeters(degToRad(90)), RADIUS_POLE, 10);
			for (int i = -90; i < 90; i += 10)
				Assert::IsTrue(earthRadiusInMeters(degToRad(i))>=RADIUS_POLE
				&& earthRadiusInMeters(degToRad(i)) <= RADIUS);
		}

		TEST_METHOD(testLocationToPoint)
		{
			Coord a = { 47.213298, -1.554762, 10.0 };
			Coord b = { 62.094063, 1.691526, 50.0 };
			Coord transfo;
			Point p;
			Vector norm;
			double radius;

			locationToPoint(&a, &p, &norm, &radius);
			transfo.elevation = a.elevation;
			pointToLocation(&transfo, &p, &radius);

			Assert::AreEqual(a.elevation, transfo.elevation, 0.001);
			Assert::AreEqual(a.latitude, transfo.latitude, 0.001);
			Assert::AreEqual(a.longitude, transfo.longitude, 0.001);

			locationToPoint(&b, &p, &norm, &radius);
			transfo.elevation = b.elevation;
			pointToLocation(&transfo, &p, &radius);

			Assert::AreEqual(b.elevation, transfo.elevation, 0.001);
			Assert::AreEqual(b.latitude, transfo.latitude, 0.001);
			Assert::AreEqual(b.longitude, transfo.longitude, 0.001);
		}

		TEST_METHOD(testrotation)
		{
			Coord a = { 47.213298, -1.554762, 10.0 };
			Coord b = { 62.094063, 1.691526, 50.0 };
			Coord transfo;
			Point p, unrotated;
			double radius, aradius = 0 , bradius = 0;
			rotateGlobe(&b, &a, bradius, &p, &radius);

			undoRotation(&a, &p, &unrotated);
			transfo.elevation = b.elevation;
			pointToLocation(&transfo, &unrotated, &radius);

			Assert::AreEqual(transfo.elevation, b.elevation, 0.5);
			Assert::AreEqual(transfo.latitude, b.latitude, 0.5);
			Assert::AreEqual(transfo.longitude+a.longitude, b.longitude, 0.5);
		}

		TEST_METHOD(testnormalizeVector)
		{
			Coord a = { 47.213298, -1.554762, 10.0 };
			Coord b = { 62.094063, 1.691526, 50.0 };


			Point p1, p2;
			Vector n1, n2, norm;
			double r1, r2, rad;
			locationToPoint(&a, &p1, &n1, &r1);
			locationToPoint(&b, &p2, &n2, &r2);


			Assert::AreEqual(normalizeVectorDiff(p2, p1,&norm, &rad),0);

			Assert::AreEqual(norm.x, -0.809, 0.001);
			Assert::AreEqual(norm.y, 0.123, 0.001);
			Assert::AreEqual(norm.z, 0.574, 0.001);
			Assert::AreEqual(rad, 1.0);

			Assert::AreEqual(normalizeVectorDiff(p1, p1, &norm, &rad), 1);
		}

		TEST_METHOD(testDistance)
		{
			Coord a = { 47.213298, -1.554762, 10.0 };
			Coord b = { 62.094063, 1.691526, 50.0 };


			Point p1, p2;
			Vector n1, n2;
			double r1, r2;
			locationToPoint(&a, &p1, &n1, &r1);
			locationToPoint(&b, &p2, &n2, &r2);

			double d = distance(p1, p2);

			Assert::AreEqual(d, 1664385, 5.0);
		}

		TEST_METHOD(testCalculate)
		{
			Coord a = { 47.213298, -1.554762, 10.0 };
			Coord b = { 48.094063, -1.691526, 50.0 };

			double distance, azimuth, altitude;
			calculate(a, b, &distance, &azimuth, &altitude);

			Assert::AreEqual(distance, 98.46, 0.01);
			Assert::AreEqual(azimuth, 354.06, 0.01);
			Assert::AreEqual(altitude, 0.42, 0.01);

		}

		TEST_METHOD(testStep)
		{
			Coord a = { 47.213298, -1.554762, 10.0 };
			Coord b = { 48.094063, -1.691526, 50.0 };

			double distance, azimuth, altitude;
			calculate(a, b, &distance, &azimuth, &altitude);

			Coord step;
			double distance_travelled = 0, tmp_dist = 0;

			double c = azimuth;
			double speed = 3.6;
			double time = 60;

			goForward(a, &step, c, speed, time);

			calculate(a, step, &tmp_dist, &c, &altitude);
			distance_travelled += tmp_dist;
			
			Assert::AreEqual(distance_travelled, speed/3600*time, 0.01);
			Assert::AreEqual(azimuth, c, 0.01);
			Assert::AreEqual(step.elevation, 10.00, 0.01);
			Assert::AreEqual(step.latitude, 47.2138, 0.01);
			Assert::AreEqual(step.longitude, -1.55484, 0.01);
		}
	};
}
