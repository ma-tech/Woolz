#include <stdio.h>
#include <Wlz.h>

int		main(int argc, char *argv[])
{
  int		squashed;
  double	dd,
  		thr;
  double	xTr[3],
  		yTr[3];
  WlzDVertex2	sVx[3],
  		dVx[3];


  thr = 1.0e-6;
  sVx[0].vtX = 260; sVx[0].vtY = 60;
  sVx[1].vtX = 240; sVx[1].vtY = 60;
  sVx[2].vtX = 250; sVx[2].vtY = 40;
  dVx[0].vtX = 19;  dVx[0].vtY = 175;
  dVx[1].vtX = 17;  dVx[1].vtY = 190;
  dVx[2].vtX = 5;   dVx[2].vtY = 178;
  (void )printf("Results should be\n");
  (void )printf("  dd = 400\n");
  (void )printf("  squashed = 0\n");
  (void )printf("  xTr = 0.1, 0.65, -46\n");
  (void )printf("  yTr = -0.75, 0.225, 356.5\n");
  (void )printf("Results are\n");
  dd = WlzGeomTriangleSnArea2(sVx[0], sVx[1], sVx[2]);
  (void )printf("  dd = %g\n", dd);
  squashed = WlzGeomTriangleAffineSolve(xTr, yTr, dd, sVx, dVx, thr);
  (void )printf("  squashed = %d\n", squashed);
  (void )printf("  xTr = %g, %g, %g\n", xTr[0], xTr[1], xTr[2]);
  (void )printf("  yTr = %g, %g, %g\n", yTr[0], yTr[1], yTr[2]);
  exit(0);
}

