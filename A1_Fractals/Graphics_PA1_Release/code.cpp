#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "code.h"
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>

//////////////////////////////////////////////////////
// The Global Variables To Be Used

struct Triangle
{
    double vertices[3][2];
    double matrix[3][3];
    int color_index;
};

// list of triangles, each element is of type Triangle
// to access the kth triangle, just use triangles[k]
// to get the number of triangles in the list, use triangles.size()
std::vector<Triangle> triangles;

// a triangle object for temporary storage of points
Triangle triangle_to_draw;
// count the number of points specified in this triangle
int point_count = 0;

// depth of recursion for IFS
// inital value is 8
// change the initial value here
int recursion_depth = 8;

// color array for triangles
// size is 11. So color_index should range from 0 to 10 for triangles.
double color_array[][3] = {
    {0.9, 0, 0}, // red
    {0, 0.5, 0.4},
    {0.1, 0.2, 0.46},
    {0.9, 0.9, 0},
    {0, 1.0, 0},
    {0, 1.0, 1.0},
    {0, 0, 1.0},
    {1.0, 0, 1.0},
    {0.9, 0.6, 0},
    {0.9, 1.0, 0.6},
    {0.2, 0.2, 0.2}};

/////////////////////////////////////////////////////////////////////////////////////////
// This function is called to clear triangles.

void ClearTriangles()
{
    triangles.clear();
    point_count = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function is called to draw triangles in the Triangles window.

void DrawTriangles()
{
    // uncomment this line if you would like to unfill triangles
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // sample code to draw vertices of triangle
    // glColor3d(1.0, 1.0, 1.0);
    // glPointSize(4);
    // glBegin(GL_POINTS);
    // for (int i = 0; i < point_count; i++)
    // {
    //     glVertex2d(triangle_to_draw.vertices[i][0], triangle_to_draw.vertices[i][1]);
    // }
    // glEnd();

    // TODO: Add code to draw triangles here. Use triangles.size() to get number of triangles.
    for (int i = 0; i < triangles.size(); i++)
    {
        glBegin(GL_TRIANGLES); // Draw filled triangle
        for (int j = 0; j < 3; j++)
        {
            // Use triangles[i] to get the ith triangle.
            // Remember to set current gl color from the color_array
            glColor3f(color_array[triangles[i].color_index][0], color_array[triangles[i].color_index][1], color_array[triangles[i].color_index][2]); // Specify color for each vertex of triangle
            glVertex3f(triangles[i].vertices[j][0], triangles[i].vertices[j][1], 0);                                                                 // Specify position for each vertex of triangle
        }
        glEnd();
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function draws the factal in a recursion manner. The depth of recursion is k.
// When k=0, the function draw the most basic shape, that is, the original triangle, and
// stops recursion.

void RecursiveFractal(int k)
{
    // TODO: Add code to implement the IFS method to draw fractal.
    // You can follow the pseudo code in the slides.
    //  for each transformation matrix Mi
    if (k > 0)
    {
        for (int i = 0; i < triangles.size(); i++)
        {
            // model view matrix = [
            //      a, d, 0, 0
            //      b, e, 0, 0
            //      0, 0, 1, 0
            //      c, f, 0, 1
            // ]

            // Mi is transformation matrix from T0 to Ti
            GLdouble temp[16] = {triangles[i].matrix[0][0], triangles[i].matrix[1][0], 0, 0,
                                 triangles[i].matrix[0][1], triangles[i].matrix[1][1], 0, 0,
                                 0, 0, 1, 0,
                                 triangles[i].matrix[0][2], triangles[i].matrix[1][2], 0, 1};
            // 1: push current modelview matrix into stack
            glPushMatrix();

            // 2: multiply current matrix with Mi
            glMultMatrixd(temp);

            // 3: Recursive Fractal (k-1) to draw a fractal in Tk-1,i
            RecursiveFractal(k - 1);
            // 4: pop matrix from matrix stack
            glPopMatrix();
        }
    }
    else
    {
        // else draw triangle T(p1,p2,p3) with current modelview matrix
        glBegin(GL_TRIANGLES);
        // glVertex2d(-0.5, -0.5);
        // glVertex2d(0.5, -0.5);
        // glVertex2d(0.0, 0.5);
        glVertex2d(triangles[0].vertices[0][0], triangles[0].vertices[0][1]);
        glVertex2d(triangles[0].vertices[1][0], triangles[0].vertices[1][1]);
        glVertex2d(triangles[0].vertices[2][0], triangles[0].vertices[2][1]);
        glEnd();
    }

    // Use triangles.size() to get number of triangles.
    // Use triangles[i] to get the ith triangle.
    // The fields of struct Triangle include:
    //	 double vertices[3][2];
    //	 double matrix[3][3];
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function invokes RecursiveFractal()

void ConstructiveFractals()
{

    if (triangles.size() < 2)
        return;

    glColor3f(1.0, 1.0, 0.0);
    RecursiveFractal(recursion_depth);
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function is called to handle mouse left click events.
// m_x and m_y are the x and y coordinates of the clicked point in OpenGL coordinate system.

void MouseInteraction(double m_x, double m_y)
{
    // TODO: Store the point temporarily into the variable triangle_to_draw.
    triangle_to_draw.vertices[point_count % 3][0] = m_x;
    triangle_to_draw.vertices[point_count % 3][1] = m_y;

    triangle_to_draw.matrix[0][point_count % 3] = 0;
    triangle_to_draw.matrix[1][point_count % 3] = 0;
    triangle_to_draw.matrix[2][point_count % 3] = 0;

    // When 3 points are specified, we get a new triangle.
    if (point_count % 3 == 2)
    {
        // if this is not the first triangle
        if (triangles.size() > 0)
        {
            // Compute the matrix for affine transformation from the first triangle to this one
            //	 by invoking AffineMatricesCalculation().
            AffineMatricesCalculation(triangles[0].vertices, triangle_to_draw.vertices, triangle_to_draw.matrix);
        }
        triangle_to_draw.color_index = triangles.size() % 11;
        // Store both the points and matrix of the triangle into a new element of the list 'triangles'
        // Use push_back function from std::Vector to add new triangles to the list
        triangles.push_back(triangle_to_draw);
    }

    point_count++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is a tool function that computes matrix for the affine transformation from one original
//   triangle to another triangle.
// Input:
//		 v_original[][2]	the pointer to an array containing data of original triangle
//		 v_transformed[][2]	the pointer to an array containing data of triangle obtained by transforming original triangle
// Output:
//		 matrix[][3]		a pointer to 3x3 matrix that is the affine transformation that changes original triangle to the later one.

void AffineMatricesCalculation(double v_original[][2], double v_transformed[][2], double matrix[][3])
{
    // TODO: Compute the affine transformation matrix that transforms triangle specified in v_original to the one specified in v_transformed.
    // If you do not want to calculate the inverse of T yourself, we provide a tool function InverseMatrix(). This function could compute the inverse of T.
    double inverse_t[3][3];

    // v_original = [
    //  [0][0]x1, [0][1]y1,
    //  [1][0]x2, [1][1]y2
    //  [2][0]x3, [2][1]y3
    // ]

    // ogMatrix = [
    //  [0][0]x1, [0][1]x2, [0][2]x3
    //  [1][0]y1, [1][1]y2, [1][2]y3
    //  [2][0]z1, [2][1]z2, [2][2]z3 = 1
    // ]
    double ogMatrix[3][3] = {
        {
            v_original[0][0],
            v_original[1][0],
            v_original[2][0],
        },
        {
            v_original[0][1],
            v_original[1][1],
            v_original[2][1],
        },
        {1, 1, 1}};
    InverseMatrix(ogMatrix, inverse_t);
    double tMatrix[3][3] = {
        {
            v_transformed[0][0],
            v_transformed[1][0],
            v_transformed[2][0],
        },
        {
            v_transformed[0][1],
            v_transformed[1][1],
            v_transformed[2][1],
        },
        {1, 1, 1}};
    // Base the computation on the formula M = T'T^(-1), where T' is the 3x3 matrix with each column the homogeneous coordinates of transformed triangle's vertex
    // and T is 3x3 matrix organized in a similar manner but stores data of the original triangle.

    // v_transformed = [
    //  [0][0]x1, [1][0]x2, [2][0]x3
    //  [0][1]y1, [1][1]y2, [2][1]y3
    //  1, 1, 1
    // ]

    // matrix pattern = [
    //  0000 + 0110 + 0220      0001 + 0111 + 0221      0002 + 0112 + 0222
    //  1000 + 1110 + 1220      1001 + 1111 + 1221      1002 + 1112 + 1222
    //  0, 0, 1
    // ]

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            matrix[i][j] = 0;

            for (int k = 0; k < 3; k++)
            {
                matrix[i][j] += tMatrix[i][k] * inverse_t[k][j];
            }
        }
    }

    // matrix[0][0] = (v_transformed[0][0] * inverse_t[0][0]) + (v_transformed[1][0] * inverse_t[1][0]) + (v_transformed[2][0] * inverse_t[2][0]);
    // matrix[0][1] = (v_transformed[0][0] * inverse_t[0][1]) + (v_transformed[1][0] * inverse_t[1][1]) + (v_transformed[2][0] * inverse_t[2][1]);
    // matrix[0][2] = (v_transformed[0][0] * inverse_t[0][2]) + (v_transformed[1][0] * inverse_t[1][2]) + (v_transformed[2][0] * inverse_t[2][2]);
    // matrix[1][0] = (v_transformed[0][1] * inverse_t[0][0]) + (v_transformed[1][1] * inverse_t[1][0]) + (v_transformed[2][1] * inverse_t[2][0]);
    // matrix[1][1] = (v_transformed[0][1] * inverse_t[0][1]) + (v_transformed[1][1] * inverse_t[1][1]) + (v_transformed[2][1] * inverse_t[2][1]);
    // matrix[1][2] = (v_transformed[0][1] * inverse_t[0][2]) + (v_transformed[1][1] * inverse_t[1][2]) + (v_transformed[2][1] * inverse_t[2][2]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// A routine to calculate inverse matrix of a 3x3 matrix which has all its values in the third row being 1.
//	original_m: 3x3 matrix with original_m[2][0]=original_m[2][1]=original_m[2][2]=1.
//  inverse_m:  3x3 matrix, the inverse of original_m.
//
void InverseMatrix(double original_m[][3], double inverse_m[][3])
{
    double determinant;
    determinant = original_m[0][0] * (original_m[1][1] - original_m[1][2]) - original_m[0][1] * (original_m[1][0] - original_m[1][2]) +
                  original_m[0][2] * (original_m[1][0] - original_m[1][1]);

    inverse_m[0][0] = (original_m[1][1] - original_m[1][2]) / determinant;
    inverse_m[1][0] = (original_m[1][2] - original_m[1][0]) / determinant;
    inverse_m[2][0] = (original_m[1][0] - original_m[1][1]) / determinant;

    inverse_m[0][1] = (original_m[0][2] - original_m[0][1]) / determinant;
    inverse_m[1][1] = (original_m[0][0] - original_m[0][2]) / determinant;
    inverse_m[2][1] = (original_m[0][1] - original_m[0][0]) / determinant;

    inverse_m[0][2] = (original_m[0][1] * original_m[1][2] - original_m[0][2] * original_m[1][1]) / determinant;
    inverse_m[1][2] = (original_m[0][2] * original_m[1][0] - original_m[0][0] * original_m[1][2]) / determinant;
    inverse_m[2][2] = (original_m[0][0] * original_m[1][1] - original_m[0][1] * original_m[1][0]) / determinant;
}
