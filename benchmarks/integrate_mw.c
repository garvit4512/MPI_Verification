/* FEVS: A Functional Equivalence Verification Suite for High-Performance
 * Scientific Computing
 *
 * Copyright (C) 2010, Stephen F. Siegel, Timothy K. Zirkel,
 * University of Delaware
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301 USA.
 */


#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<math.h>
#include<float.h>
#include<mpi.h>

#define WORK_TAG 1
#define TERMINATE_TAG 2
#define RESULT_TAG 3

#pragma TASS input int
#define INTERVALS 300
#pragma TASS input
double a;
#pragma TASS input
double b;
#pragma TASS input
double tolerance;

double epsilon = DBL_MIN*10;
double pi = M_PI;
int nprocs, rank, numWorkers, numTasks;
double subintervalLength, subtolerance;

#pragma TASS output
double out1;

#pragma TASS abstract continuous(100) double sin(double x);
//#pragma TASS abstract continuous(100) double cos(double x);

#pragma TASS abstract continuous(0) double fabs(double x);
/*
double fabs(double x) {
	if (x >= 0.0) {
		return x;
	} else {
		return -x;
	}
}*/

//double f(double x) { return sin(x); }
#pragma TASS abstract continuous(0) double f(double x);

//double f2(double x) { return cos(x); }

/* Sequential function to recursively estimate integral of f from a to
 * b.  fa=f(a), fb=f(b), area=given estimated integral of f from a to
 * b.  The interval [a,b] is divided in half and the areas of the two
 * pieces are estimated using Simpson's rule.  If the sum of those two
 * areas is within tolerance of given area, convergence has been
 * achieved and the result returned is the sum of the two areas.
 * Otherwise the function is called recursively on the two
 * subintervals and the sum of the results returned.
 */
double integrate_seq(double a, double b, double fa, double fb, double area,
		     double tolerance) {
  double delta = b - a;
  double c = a+delta/2;
  double leftArea, rightArea, fc;
  double absolute, leftIntegral, rightIntegral;
  fc = f(c);
  leftArea = (fa+fc)*delta/4;
  rightArea = (fc+fb)*delta/4;

  if (tolerance < epsilon) {
    //printf("Tolerance may not be possible to obtain.\n");
    return leftArea+rightArea;
  }
  absolute = fabs(leftArea+rightArea-area);
  if (absolute<=tolerance) {
    return leftArea+rightArea;
  }
  leftIntegral = integrate_seq(a, c, fa, fc, leftArea, tolerance/2);
  rightIntegral = integrate_seq(c, b, fc, fb, rightArea, tolerance/2);
  return leftIntegral + rightIntegral;
}

/* Sequential algorithm to estimate integral of f from a to b within
 * given tolerance. */
double integral_seq(double a, double b, double tolerance) {
  double fa, fb, result;
  fa = f(a);
  fb = f(b);
  result = integrate_seq(a, b, fa, fb, (fa+fb)*(b-a)/2, tolerance);
  return result;
}

/* Worker function: called by procs of non-0 rank only.  Accepts
 * task from manager.   Each task is represented by a single integer
 * i, corresponding to the i-th subinterval of [a,b].  The task is
 * to compute the integral of f(x) over that subinterval within
 * tolerance of subtolerance.   At the end, a Bcast is used to
 * get final answer from manager. */
double worker(double a) {
  double left, right, result, answer;
  int task;
  int control = 1;
  MPI_Status status;

  // while (control != 0) {
    MPI_Recv(&task, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    //if (status.MPI_TAG == WORK_TAG) {
      left = a+task*subintervalLength;
      right = left+subintervalLength;
      result = integral_seq(left, right, subtolerance);
      MPI_Send(&result, 1, MPI_DOUBLE, 0, RESULT_TAG, MPI_COMM_WORLD);
      //} else {
      //control = 0;
      // }
      //}
  //fflush(stdout);
  return 0.0;
}

/* manager: this function called by proc 0 only.  Distributes tasks
 * to workers and accumulates results.   Each task is represented
 * by an integer i.  The integer i represents a subinterval
 * of the interval [a,b].   Returns the final result.  */
double manager() {
  int i, task = 0;
  MPI_Status status;
  double result, answer = 0.0;
    
  for (i=1; i<nprocs; i++) {
    MPI_Send(&task, 1, MPI_INT, i, WORK_TAG, MPI_COMM_WORLD);
    task++;
  }
  while (task < numTasks) {
    MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, RESULT_TAG,
	     MPI_COMM_WORLD, &status);
    //MPI_Send(&task, 1, MPI_INT, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD);
    //  task++;
    // answer += result;
    for(i=1;i<nprocs; i++){
      if (task < numTasks){
	MPI_Send(&task, 1, MPI_INT, i, WORK_TAG, MPI_COMM_WORLD);
	task++;
      }
      //if(task == numTasks) break;
    }

  }
  /* for (i=1; i<nprocs; i++) { */
  /*   MPI_Recv(&result, 1, MPI_DOUBLE, i, RESULT_TAG, */
  /* 	     MPI_COMM_WORLD, &status); */
  /*   MPI_Send(NULL, 0, MPI_INT, status.MPI_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD); */
  /*   answer += result; */
  /* } */

  for (i=1; i<=numTasks; i++) {
    MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, RESULT_TAG,
	     MPI_COMM_WORLD, &status);
    /* MPI_Send(NULL, 0, MPI_INT, status.MPI_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD); */
    answer += result;
  }

  return answer;
}

/* called in collective fashion */
double integral(double a, double b, double tolerance) {
  double result;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  numWorkers = nprocs - 1;
  //  numTasks = INTERVALS;
  numTasks = numWorkers;
//  assert (numTasks == numIntervals);
  subintervalLength = (b-a)/numTasks;
  subtolerance = tolerance/numTasks;
  if (rank == 0) {
    result = manager();
    return result;
  } else {
    result = worker(a);
    return result;
  }
}

int main() {
  int rank, argc;
  char **argv;
  double localResult;
#pragma TASS assume b-a <= 1.0 && b-a >= 0.0;
#pragma TASS assume forall{double x} fabs(x) == x || fabs(x) == -x;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  localResult = integral(a, b, tolerance);
  if (rank == 0)
    out1 = localResult;
  //if (rank == 0) { printf("%4.20lf\n", result); }
  //result = integral(0, pi, .0000000001, f2);
  //if (rank == 0) { printf("%4.20lf\n", result); }
  MPI_Finalize();

  return 0;
}
