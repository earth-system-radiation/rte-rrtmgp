
#pragma once

#include "const.h"


inline void reorder_123x321_kernel(int d1, int d2, int d3, real3d const &array_in, real3d &array_out) {
  int constexpr TILE_SIZE = 2;
  int ntiles1 = d1 / TILE_SIZE + 1;
  int ntiles3 = d3 / TILE_SIZE + 1;

  // for (int i2=1; i2<=d2; i2++) {
  //   for (int t3=1; t3<=ntiles3; t3++) {
  //     for (int t1=1; t1<=ntiles1; t1++) {
  //       for (int it3=1; it3<=TILE_SIZE; it3++) {
  //         for (int it1=1; it1<=TILE_SIZE; it1++) {
  parallel_for( Bounds<5>(d2,ntiles3,ntiles1,TILE_SIZE,TILE_SIZE) , YAKL_LAMBDA (int i2, int t3, int t1, int it3, int it1) {
    int i3 = (t3-1)*TILE_SIZE + it3;
    int i1 = (t1-1)*TILE_SIZE + it1;
    if (i3 <= d3 && i1 <= d1) {
      array_out(i3,i2,i1) = array_in(i1,i2,i3);
    }
  });
}


inline void reorder_123x312_kernel(int d1, int d2, int d3, real3d const &array_in, real3d &array_out) {
  // do i3 = 1 , d3
  //   do i2 = 1 , d2
  //     do i1 = 1 , d1
  parallel_for( Bounds<3>(d3,d2,d1) , YAKL_LAMBDA (int i3, int i2, int i1) {
    array_out(i3,i1,i2) = array_in(i1,i2,i3);
  });
}


