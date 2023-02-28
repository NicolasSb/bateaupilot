#include "command.h"

void setupCommand()
{
	//for simulation purpose
	srand((unsigned int)time(NULL));   // Initialization, should only be called once.

	//init serial com

	//configure input output

	//setup interrupt routines

	//write default motor command

	//another interrupt to compute the PID correction (-> servoing)

}

float command()
{
	return 0.0;
}


void perturbation(Coord *c)
{
	//drift effect / gps imprecision, ...

	int exp = rand() % 2;
	double delta = (float)pow(-1, exp) * (float)rand() / (float)RAND_MAX;
	double eta = (float)pow(-1, exp) * (float)rand() / (float)RAND_MAX;
	c->latitude += delta/10000;
	c->longitude += eta/10000;	
}

void servoing(double target, double measure, double *command)
{
	t += SAMPLING_FREQ;

	error = target - measure;

	errorsum += error;

	int vmot = (int)(KP * error + KD * (error - last_error) + KI * errorsum);
	//normalize and control motor
	//if (vmot > 255) vmot = 255;
	//else if (vmot < -255) vmot = -255;

	*command = (double)vmot + target ;

	//command(vmot);
}

float normalizeAngle(float angle)
{
	float ret = angle;
	while (ret < 0 || ret > 360)
	{
		if (ret > 360) ret -= 360;
		if (ret < 0) ret += 360;
	}
	return ret;
}