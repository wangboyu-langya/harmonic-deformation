#include "Reduction.h"
#include <igl/convex_hull.h>

using namespace igl;
using namespace Eigen;
using namespace std;

void Reduction(Eigen::MatrixXd &V, Eigen::MatrixXd &F, int dim) {
    // reduce one dimension of V
    VectorXd x;
    x.setZero(V.rows());
    V.col(dim) = x;
    // find the boundary of reduced point cloud and triangulate it.
    convex_hull(V, F);
}
