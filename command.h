#ifndef COMMAND_H
#define COMMAND_H

#include <time.h>
#include <stdlib.h>
#include "coord.h"

#define KP 0.9
#define KI 0
#define KD 0

#define SAMPLING_FREQ 20
/* define input output ports*/
//const int pinOutput1 = x;		//

/*define sampling variables*/
static unsigned int t = 0;

static double error = 0.0;
static double last_error = 0.0;
static double errorsum = 0.0;

void setupCommand();

void stepMotor(int i);

double commandMotor(double target);

void loop();

void servoing(double target, double measure, double * command);

float command();

void perturbation(Coord * c);

float normalizeAngle(float angle);

#endif // !COMMAND_H
