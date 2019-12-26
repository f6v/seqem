#ifndef __MAX_H__
#define __MAX_H__

const int max_rdln = 2048;
const int max_name = 80;
const int max_iteration_offset_for_matrix = 5;							// otherwise we get a zero denominator which isn't good
const int min_iter = 2 * max_iteration_offset_for_matrix;
//const int max_iter = 200 + min_iter;
const int max_matrix_iter = 20;
const int min_slice = 20;
#endif
