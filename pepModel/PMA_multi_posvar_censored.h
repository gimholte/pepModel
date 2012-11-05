/*
 * PMA_multi_posvar_censored.h
 *
 *  Created on: Nov 2, 2012
 *      Author: Gregory Imholte
 */

#ifndef PMA_MULTI_POSVAR_CENSORED_H_
#define PMA_MULTI_POSVAR_CENSORED_H_

struct MH_TUNE{
	int count;
	int total_count;
	double tune;
};

typedef struct MH_TUNE adpt;

#endif /* PMA_MULTI_POSVAR_CENSORED_H_ */
