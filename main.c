#include <stdlib.h>
#include <stdio.h>
#include <windows.h>

#include "coord.h"
#include "command.h"

int main()
{
	setupCommand();
	Coord a = { 47.213298, -1.554762, 10.0 };
	Coord b = { 48.094063, -1.691526, 50.0 };

	printf("Coordonnees de chez Blandine %lf, %lf\n", a.latitude, a.longitude);
	printf("Coordonnees de chez Nico %lf, %lf\n", b.latitude, b.longitude);

	double distance, azimuth, altitude;
	calculate(a, b, &distance, &azimuth, &altitude);

	printf("distance : % lf\n", distance);
	printf("Azimuth : % lf\n", azimuth);
	printf("Altitude : % lf \n", altitude);

	Coord step;	
	double distance_travelled = 0, tmp_dist = 0;

	double c = azimuth;
	double measure;
	//double diff = 0;
	//correcteur ?
	//actionneur ? 
	//perturbation ? 
	//sortie
	double speed = 36;
	double time = 60;

	printf("\n\n");

	while (distance - distance_travelled > 0)
	{
		goForward(a, &step, c, speed, time);
		perturbation(&step);

		calculate(a, step, &tmp_dist, &measure, &altitude);
		distance_travelled += tmp_dist;
		a = step;


		calculate(a, b, &distance, &azimuth, &altitude);

		servoing(azimuth, measure, &c);

		//diff = perturbation();

		printf("distance : % lf\n", distance_travelled);
		printf("Azimuth : %f - target : %lf\n", c, azimuth);
		Sleep(1000);
	}

	return 0;
}