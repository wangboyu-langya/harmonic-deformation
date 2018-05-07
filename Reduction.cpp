#include <igl/triangle/triangulate.h>
#include <igl/boundary_facets.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unique.h>
#include "convex_hull.cpp"
using namespace igl;
using namespace Eigen;
using namespace std;

void Reduction(Eigen::MatrixXd &V, Eigen::MatrixXd &F, int dim) {
    // parameter of cage
    double factor = 1.2;
    // reduce one dimension of V
    VectorXd x;
    x.setZero(V.rows());
    // reduce V to 2d
    if (dim == 0) {V.col(0) = V.col(2); V.col(2) = x;}
    if (dim == 1) {V.col(1) = V.col(2); V.col(2) = x;}
	NoChange_t change = NoChange;
    V.conservativeResize(change, 2);
    // find the inner boundary of reduced point cloud
    Eigen::MatrixXd EI, BI;
    convex_hull(V, EI, BI);
    // get the outer boundary
    Eigen::MatrixXd BO, EO, VB;
    BO.resize(BI.rows(), 2);
    for (int i = 0; i < BO.rows(); ++i) {
        for (int j = 0; j < BO.cols(); ++j) {
            if (BI(i, j) > 0) BO(i, j) = BI(i, j) * factor;
            else BO(i, j) *= - BI(i, j) * factor;
        }
    }
    EO = EI.replicate(1, 1);
    // get the outer triangularization
    Eigen::MatrixXd E, B;
    E.resize(2 * EI.rows(), 2);
    B.resize(2 * BI.rows(), 2);
    E << EI, EO;
    B << BI, BO;
    Eigen::MatrixXd VO, FO;
    Eigen::RowVector2d y;
    y << 0, 0;
    for (int i = 0; i < BI.rows(); ++i) {
        y += BI.row(i);
    }
    y(0) /= BI.rows();
    y(1) /= BI.rows();
    igl::triangle::triangulate(B, E , y, "a0.005q", VO, FO);
    // get the inner triangularization
    Eigen::MatrixXd VI, FI, HI;
    igl::triangle::triangulate(BI, EI , HI, "a0.005q", VI, FI);
    // get the final V
    MatrixXd VF, FF;
    VF.resize(VO.rows() + VI.rows(), 2);
    FF.resize(FO.rows() + FI.rows(), 2);
    VF << VI, VO;
    FF << FI, FO;
//    BI.conservativeResize(2 * BI.rows(), NoChange_t);
//    Eigen::VectorXi RB;
//    RB.setLinSpaced(BI.rows(), 0 + BI.rows(), BI.rows() * 2 - 1);
//    igl::slice_into(BO, RB, 1, BI);
    // Triangulate the interior
    // Find boundary vertex
    Eigen::MatrixXi EF;
    Eigen::VectorXi b, IA, IC;
    igl::boundary_facets(FF, EF);
    igl::unique(EF, b, IA, IC);
    //cage vertex array, #V x3
    Eigen::MatrixXd VC(0, 2);
    igl::slice(VF, b, 1, VC);
    //get boundary edge
    VectorXd R1, R2;
    R1.setLinSpaced(VC.rows(), 0, VC.rows() - 1);
    R2.setLinSpaced(VC.rows(), 1, VC.rows());
    R2(VC.rows() - 1) = 0;
    Eigen::MatrixXd C1, C2;
    slice(VC, R1, 1, C1);
    slice(VC, R2, 1, C2);
    // plot the cage
//    igl::opengl::glfw::Viewer viewer;
//    viewer.data().clear();
//    viewer.data().set_mesh(VF, FI);
//    viewer.data().point_size = 20;
//    viewer.data().add_points(VC, Eigen::RowVector3d(1, 0, 0));
//    viewer.data().add_edges(C1, C2, Eigen::RowVector3d(1, 0, 0));
}
