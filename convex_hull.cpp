// Inputs:
//   V #V by 2 list of 2D vertex positions
// Outputs:
//   E #E by 2 list of vertex ids forming unoriented edges of the boundary of the polygon
// #include "convex_hull.h"
#include <vector>
#include <Eigen/Dense>

using namespace std;
//int orientation(Eigen::MatrixBase<DerivedV> & p, Eigen::MatrixBase<DerivedE> & q, Eigen::MatrixBase<DerivedE> & r)
int orientation(Eigen::RowVector2d &p, Eigen::RowVector2d &q, Eigen::RowVector2d &r)
{
	double val = (q(1) - p(1)) * (r(0) - q(0)) - (q(0) - p(0)) * (r(1) - q(1));
	return (val > 0)? 1: 2; // clock or counterclock wise
}

void convex_hull(Eigen::MatrixXd & V, Eigen::MatrixXd & E, Eigen::MatrixXd & B){
    // Initialize Result
    vector<int> hull;
    int n = V.rows();

    // Find the leftmost point
    int l = 0;
    for (int i = 1; i < V.rows(); i++)
        if (V(i, 0) < V(l, 0))
            l = i;

    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    int p = l, q;
    do
    {
        // Add current point to result
        hull.push_back(p);
        // Search for a point 'q' such that orientation(p, x,
        // q) is counterclockwise for all points 'x'. The idea
        // is to keep track of last visited most counterclock-
        // wise point in q. If any point 'i' is more counterclock-
        // wise than q, then update q.
        Eigen::RowVector2d pv = V.row(p);
        q = (p+1)%n;
        Eigen::RowVector2d qv = V.row(q);
        q = (p+1)%n;
        for (int i = 0; i < n; i++)
        {
            // If i is more counterclockwise than current q, then
            // update q
            Eigen::RowVector2d iv = V.row(i);
//            if (orientation(V.row(p), V.row(i), V.row(q)) == 2)
            if (orientation(pv, iv, qv) == 2)
                q = i;
        }
        // Now q is the most counterclockwise with respect to p
        // Set p as q for next iteration, so that q is added to
        // result 'hull'
        p = q;
    } while (p != l);  // While we don't come to first point
    // get the boundary
    B.resize(hull.size(), 2);
    for (int i = 0; i < hull.size(); ++i)
        B.row(i) = V.row(hull[i]);
    // get the edge
    E.resize(hull.size() - 1, 2);
    for (int i = 0; i < hull.size() - 2; ++i) {
        E(i, 0) = hull[i];
        E(i, 1) = hull[i + 1];
    }
}

