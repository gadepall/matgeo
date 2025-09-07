/* points.c */
#include <math.h>
#include <stdio.h>

/* Compute shortest distance between two lines:
   Line1: point P1=(0,0,1), direction d1=(1,1,-1)
   Line2: point P2=(0,0,0), direction d2=(1,0,1)
*/
void compute_distance(double *result) {
  /* Points */
  double p1x = 0.0, p1y = 0.0, p1z = 1.0;
  double p2x = 0.0, p2y = 0.0, p2z = 0.0;

  /* Directions */
  double d1x = 1.0, d1y = 1.0, d1z = -1.0;
  double d2x = 1.0, d2y = 0.0, d2z = 1.0;

  /* P2 - P1 */
  double rx = p2x - p1x;
  double ry = p2y - p1y;
  double rz = p2z - p1z;

  /* cross = d1 x d2 */
  double cx = d1y * d2z - d1z * d2y;
  double cy = d1z * d2x - d1x * d2z;
  double cz = d1x * d2y - d1y * d2x;

  double numer = fabs(rx * cx + ry * cy + rz * cz);
  double denom = sqrt(cx * cx + cy * cy + cz * cz);

  *result = numer / denom;
}
