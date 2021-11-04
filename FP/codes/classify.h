#pragma once
#ifndef CLASSIFY_H
#define CLASSIFY_H

#include "common.h"

Str dist(Str Out, Ma_Str_M ref, Ma_Str_M que, Str GenoID, Ma_Str_X taxo, double tthr, double sthr, Ma_Str_Pa_Str_Str check);

double SD_distance(Str SD, Ve_Ma_C_D que, Ma_C_D pb, int len);

double TA_distance(Ve_Ma_C_D ref, Ma_C_D rpb, Ve_Ma_C_D que, Ma_C_D qpb, int len);

#endif