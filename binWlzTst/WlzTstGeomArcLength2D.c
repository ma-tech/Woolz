#include <stdio.h>
#include <float.h>
#include <Wlz.h>

int		main(int argc, char *argv[])
{
  double	d;
  WlzDVertex2	a,
  		b,
		c;

  c.vtX = 200.0;
  c.vtY = 200.0;
  a.vtX = 100.0;
  a.vtY = 200.0;
  b.vtX = 300.0;
  b.vtY = 200.0;
  d = WlzGeomArcLength2D(a, b, c); 
  (void )printf("% 4.4lf\n", d);
}
