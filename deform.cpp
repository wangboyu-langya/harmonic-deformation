#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <Eigen/Core>
#include <igl/boundary_facets.h>
#include <igl/colon.h>
#include <igl/cotmatrix.h>
#include <igl/jet.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/readOFF.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/unique.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Sparse>
#include <iostream>

using namespace std;
//vertex array, #V x3
Eigen::MatrixXd V(0, 3);
//cage vertex array, #V x3
Eigen::MatrixXd VH(0, 3);
//cage edge
Eigen::MatrixXd C1, C2;
//face array, #F x3
Eigen::MatrixXi F(0, 3);
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0, 1);
//handle vertex, #H x1
Eigen::VectorXi handle_V(0, 1);
bool ifdebug = true;

enum MouseMode {
    SELECT, TRANSLATE, ROTATE, NONE
};
MouseMode mouse_mode = NONE;

//rotation and translation for the handle being moved
Eigen::Vector3f translation(0, 0, 0);

bool callback_mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
    if (button == (int) igl::opengl::glfw::Viewer::MouseButton::Right)
        return false;

    down_mouse_x = viewer.current_mouse_x;
    down_mouse_y = viewer.current_mouse_y;

    if (mouse_mode == TRANSLATE) {
        int vi = pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
        if (vi >= 0 && handle_id[vi] >= 0)  //if a region was found, mark it for translation/rotation
        {
            moving_handle = handle_id[vi];
            get_new_handle_locations();
            doit = true;
        }
    }
    return doit;
}
int pickVertex(viewer.current_mouse_x, viewer.current_mouse_y) {
    int nearest = -1;
    double d = 0.0;
    for (int i = 0; i < VH.rows(); ++i) {
        double d =
    }
}

bool callback_mouse_move(igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y) {
    if (!doit)
        return false;
    if (mouse_mode == SELECT) {
        lasso->strokeAdd(mouse_x, mouse_y);
        return true;
    }
    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) {
        if (mouse_mode == TRANSLATE) {
            translation = computeTranslation(viewer,
                                             mouse_x,
                                             down_mouse_x,
                                             mouse_y,
                                             down_mouse_y,
                                             handle_centroids.row(moving_handle));
        } else {
            rotation = computeRotation(viewer,
                                       mouse_x,
                                       down_mouse_x,
                                       mouse_y,
                                       down_mouse_y,
                                       handle_centroids.row(moving_handle));
        }
        get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
        solve(viewer);
        down_mouse_x = mouse_x;
        down_mouse_y = mouse_y;
#endif
        return true;

    }
    return false;
}

bool callback_mouse_up(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
    if (!doit)
        return false;
    doit = false;
    if (mouse_mode == SELECT) {
        selected_v.resize(0, 1);
        lasso->strokeFinish(selected_v);
        return true;
    }

    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) {
#ifdef UPDATE_ONLY_ON_UP
        if(moving_handle>=0)
          solve(viewer);
#endif
        translation.setZero();
        rotation.setZero();
        rotation[3] = 1.;
        moving_handle = -1;

        compute_handle_centroids();

        return true;
    }

    return false;
};


int main() {

    using namespace Eigen;
    using namespace igl;
    // get the cage
    igl::readOFF("../data/cube_2d_1_cage.off", VH, F);
    assert(VH.rows() > 0);
    VectorXd R1, R2;
    R1.setLinSpaced(VH.rows(), 0, VH.rows() - 1);
    R2.setLinSpaced(VH.rows(), 1, VH.rows());
    R2(VH.rows() - 1) = 0;
    slice(VH, R1, 1, C1);
    if (ifdebug) cout << "C1: " << C1 << endl;
    slice(VH, R2, 1, C2);
    if (ifdebug) cout << "C2: " << C2 << endl;
    // load the mesh
    igl::readOFF("../data/cube_2d_1.off", V, F);
    assert(V.rows() > 0);
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    // plot the cage
    viewer.data().point_size = 20;
    viewer.data().add_points(VH, Eigen::RowVector3d(1, 0, 0));
    viewer.data().add_edges(C1, C2, Eigen::RowVector3d(1, 0, 0));

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    viewer.callback_key_down = callback_key_down;
    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Deformation Controls", ImGuiTreeNodeFlags_DefaultOpen)) {

            // Expose an enumeration type
            ImGui::Combo("MouseMode", (int *) (&mouse_mode), "TRANSLATE\0NONE\0\0");
        }
    }
        viewer.launch();
}



////function declarations (see below for implementation)
//bool solve(igl::opengl::glfw::Viewer &viewer);
//
//Eigen::MatrixXi E;
//Eigen::VectorXi b,IA,IC;
//Eigen::VectorXi all,in;
//Eigen::SparseMatrix<double> L,L_in_in,L_in_b;
//void solve_scalar(int dim, Eigen::VectorXd & Z) {
//    using namespace Eigen;
//  // Dirichlet boundary conditions from z-coordinate
//  VectorXd bc;
//  Z = V.col(dim);
//  igl::slice(Z,b,bc);
//
//  // Solve PDE
//  SimplicialLLT<SparseMatrix<double > > solver(-L_in_in);
//  VectorXd Z_in = solver.solve(L_in_b*bc);
//  // slice into solution
//  igl::slice_into(Z_in,in,Z);
//}

