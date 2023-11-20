// Raphaele Michelle Guillemot 3033093337
// COMP3271
// Programming Assignment #2

#include "MainFrame.h"
#include <iostream>
#include <math.h>
namespace
{

    float scale = 1.f;
    float aspect = 1.f;

#ifdef __APPLE__
    unsigned int SCR_WIDTH = 600;
    unsigned int SCR_HEIGHT = 600;
#else
    unsigned int SCR_WIDTH = 1000;
    unsigned int SCR_HEIGHT = 1000;
#endif

    void ScrollCallback(GLFWwindow *window, double xoffset, double yoffset)
    {
        scale *= std::pow(1.1f, (float)yoffset);
    }

    void FrameBufferSizeCallback(GLFWwindow *window, int width, int height)
    {
        SCR_WIDTH = width;
        SCR_HEIGHT = height;
        glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);       // Set the viewport to cover the new window
        aspect = (float)SCR_WIDTH / (float)SCR_HEIGHT; // Set the aspect ratio of the clipping area to match the viewport
    }

}

void MainFrame::LeftMouseMove(float start_x, float start_y, float curr_x, float curr_y)
{
    if (modeling_state_ == OBJ_ROTATION)
    {
        // ---------------------------------- Object Rotation ---------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "transform_mat" to rotate the object about its center.

        glm::mat4x4 transform_mat(1.f);

        // Step 1 : Find the rotation axis and determine the rotation angle
        glm::vec2 V = glm::vec2((curr_x - start_x), (curr_y - start_y));
        glm::vec2 A = glm::vec2((-V[1]), (V[0]));

        glm::vec3 Sstart = Screen2World(start_x, start_y, -1.f);
        glm::vec3 Sstart_A = Screen2World(start_x + A[0], start_y + A[1], -1.f);

        glm::vec3 ra = normalize(Sstart_A - Sstart);

        float angle = 0.01 * sqrt((A[0] * A[0]) + (A[1] * A[1]));
        // Step 2 : Get the rotation matrix about the rotation axis
        glm::mat4x4 R;
        R = glm::rotate(glm::mat4x4(1.f), angle, ra);

        //(a) Translate the rotation center to the origin
        transform_mat = translate(transform_mat, mesh_.center_);
        // std::cout << "To origin?\n";
        // std::cout << transform_mat[0][0] << " " << transform_mat[0][1] << " " << transform_mat[0][2] << " " << transform_mat[0][3] << "\n"
        //           << transform_mat[1][0] << " " << transform_mat[1][1] << " " << transform_mat[1][2] << " " << transform_mat[1][3] << "\n"
        //           << transform_mat[2][0] << " " << transform_mat[2][1] << " " << transform_mat[2][2] << " " << transform_mat[2][3] << "\n"
        //           << transform_mat[3][0] << " " << transform_mat[3][1] << " " << transform_mat[3][2] << " " << transform_mat[3][3] << "\n";

        //(b) Rotate round the origin
        transform_mat = transform_mat * R;
        // std::cout << "To R?\n";
        // std::cout << transform_mat[0][0] << " " << transform_mat[0][1] << " " << transform_mat[0][2] << " " << transform_mat[0][3] << "\n"
        //           << transform_mat[1][0] << " " << transform_mat[1][1] << " " << transform_mat[1][2] << " " << transform_mat[1][3] << "\n"
        //           << transform_mat[2][0] << " " << transform_mat[2][1] << " " << transform_mat[2][2] << " " << transform_mat[2][3] << "\n"
        //           << transform_mat[3][0] << " " << transform_mat[3][1] << " " << transform_mat[3][2] << " " << transform_mat[3][3] << "\n";

        // (c) Translate the rotation center back
        transform_mat = translate(transform_mat, -mesh_.center_);
        // std::cout << "To og origin?\n";
        // std::cout << transform_mat[0][0] << " " << transform_mat[0][1] << " " << transform_mat[0][2] << " " << transform_mat[0][3] << "\n"
        //           << transform_mat[1][0] << " " << transform_mat[1][1] << " " << transform_mat[1][2] << " " << transform_mat[1][3] << "\n"
        //           << transform_mat[2][0] << " " << transform_mat[2][1] << " " << transform_mat[2][2] << " " << transform_mat[2][3] << "\n"
        //           << transform_mat[3][0] << " " << transform_mat[3][1] << " " << transform_mat[3][2] << " " << transform_mat[3][3] << "\n";

        // Step 4 : Apply the transformation
        mesh_.ApplyTransform(transform_mat);
    }
    else if (modeling_state_ == OBJ_TRANSLATION)
    {
        // ---------------------------------- Object Translation ------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "trans_mat" to translate the object along the view plane.

        glm::mat4x4 transform_mat(1.f);

        // Step 1: For the start mouse position Sstart, find Pstart in world coordinate
        glm::vec3 Pstart = Screen2World(start_x, start_y, -1.f);
        // Step 2: For the current mouse position Scurr, find Pcurr in world coordinate
        glm::vec3 Pcurr = Screen2World(curr_x, curr_y, -1.f);

        // Step 3: Create the translation matrix (Pstart -> Pcurr)
        transform_mat = glm::translate(transform_mat, Pcurr - Pstart);

        // Step 4: Apply the transformation
        mesh_.ApplyTransform(transform_mat);
    }
    else if (modeling_state_ == OBJ_EXTRUDE)
    {
        // ---------------------------------- Face Extrusion ------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "trans_mat" to translate the face vertices along the face normal.

        // Step 1: For the current intersected face by ray R_0, calculate the normal N
        glm::vec3 Sstart = Screen2World(start_x, start_y);  // Sstart
        std::tuple<glm::vec3, glm::vec3> R_0 = Screen2WorldRay(start_x, start_y);   // R_0
        std::tuple<int, glm::vec3> Pstart = mesh_.FaceIntersection(std::get<0>(R_0), std::get<1>(R_0));
        int face_index = std::get<0>(Pstart);               // Face index
        glm::vec3 face_point = std::get<1>(Pstart);         // Face Point

        // get N by normalise (cross (some points on the face))
        Face face_o = mesh_.faces_[face_index];
        glm::vec3 N = glm::normalize(glm::cross(mesh_.vertices_[face_o[1]] - mesh_.vertices_[face_o[0]],
                                                mesh_.vertices_[face_o[2]] - mesh_.vertices_[face_o[1]]));

        // Step 2: Calculate Pcurr
        std::tuple<glm::vec3, glm::vec3> R_1 = Screen2WorldRay(curr_x, curr_y);     // R_1
        glm::vec3 VR_1 = std::get<1>(R_1);
        std::tuple<glm::vec3, glm::vec3> M(Sstart, N);  // M

        glm::vec3 V = glm::normalize(cross(N, VR_1));   // V
        glm::vec3 Q = glm::normalize(cross(V, VR_1));   // Plane Q

        float t = glm::dot(Q, Sstart) / glm::dot(Q, N); // Pcurr's t
        glm::vec3 Pcurr = Sstart + (t * N);             // Pcurr

        // Step 3: Create the translation matrix (Pstart -> Pcurr)
        glm::mat4x4 transform_mat(1.f);
        transform_mat = glm::translate(transform_mat, Pcurr - Sstart);

        // Step 4: Apply the transformation
        mesh_.ApplyFaceTransform(face_index, transform_mat);
    }
}

void MainFrame::VisualizeWorldSpace()
{
    // ---------------------------------- World Space Visualization ------------------------------------
    // TODO: Add your code here to visualize the world space.

    glBegin(GL_LINES);
    for (int i = -10; i <= 10; i++)
    {
        glColor3f(.7, .7, 1);
        glVertex3f(-10, i, -1);
        glVertex3f(10, i, -1);
        glColor3f(.7, .7, 1);
        glVertex3f(i, -10, -1);
        glVertex3f(i, 10, -1);
    };
    glEnd();
}

// -------------------------------------------------------------------------------------
// -------------------------- No need to change ----------------------------------------
// -------------------------------------------------------------------------------------
void MainFrame::MainLoop()
{
    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_SAMPLES, 4);

    // glfw window creation, set viewport with width=1000 and height=1000
    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "3DModeling", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, FrameBufferSizeCallback);
    glfwSetScrollCallback(window, ScrollCallback);
    // glad: load all OpenGL function pointers
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    const float alpha = 0.3f;
    const float beta = 0.1f;

    const float r = 5.f;
    camera_.LookAt(r * glm::vec3(std::cos(alpha) * std::cos(beta), std::cos(alpha) * std::sin(beta), std::sin(alpha)),
                   glm::vec3(0.f, 0.f, 0.f),
                   glm::vec3(0.f, 0.f, 1.f));

    glEnable(GL_DEPTH_TEST);

    // render loop
    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        // glEnable(GL_DEPTH_TEST);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Apply camera projection;
        camera_.Perspective(90.f, aspect, .5f, 10.f);
        camera_.UpdateScale(scale);
        scale = 1.f;
        camera_.ApplyProjection();

        glClearColor(0.0, 0.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear the display

        DrawScene();

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        glfwPollEvents();
        glfwSwapBuffers(window);
    }

    // glfw: terminate, clearing addl previously allocated GLFW resources.
    glfwTerminate();
}

void MainFrame::ProcessInput(GLFWwindow *window)
{
    // Key events
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_ROTATION;
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_TRANSLATION;
    }
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_SUBDIVIDE;
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_EXTRUDE;
    }

    int current_l_mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

    // Handle left mouse
    if (current_l_mouse_state == GLFW_PRESS)
    {
        double xposd, yposd;
        float xpos, ypos;
        glfwGetCursorPos(window, &xposd, &yposd);
        xpos = float(xposd);
        ypos = float(SCR_HEIGHT - yposd);
        if (l_mouse_state_ == GLFW_RELEASE)
        {
            LeftMouseClick(xpos, ypos);
            l_click_cursor_x_ = xpos;
            l_click_cursor_y_ = ypos;
        }
        if (l_mouse_state_ == GLFW_PRESS &&
            (std::abs(xpos - last_cursor_x_) > 2.f || std::abs(ypos - last_cursor_y_) > 2.f))
        {
            LeftMouseMove(l_click_cursor_x_, l_click_cursor_y_, xpos, ypos);
        }
        last_cursor_x_ = float(xpos);
        last_cursor_y_ = float(ypos);
    }
    if (current_l_mouse_state == GLFW_RELEASE)
    {
        if (l_mouse_state_ == GLFW_PRESS)
        {
            LeftMouseRelease();
        }
    }
    l_mouse_state_ = current_l_mouse_state;

    // Handle right mouse
    int current_r_mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
    if (current_r_mouse_state == GLFW_PRESS)
    {
        double xposd, yposd;
        float xpos, ypos;
        glfwGetCursorPos(window, &xposd, &yposd);
        xpos = float(xposd);
        ypos = float(SCR_HEIGHT - yposd);
        if (r_mouse_state_ == GLFW_RELEASE)
        {
            RightMouseClick(xpos, ypos);
        }
        if (r_mouse_state_ == GLFW_PRESS &&
            (std::abs(xpos - last_cursor_x_) > 2.f || std::abs(ypos - last_cursor_y_) > 2.f))
        {
            RightMouseMove(last_cursor_x_, last_cursor_y_, xpos, ypos);
        }
        last_cursor_x_ = float(xpos);
        last_cursor_y_ = float(ypos);
    }
    if (current_r_mouse_state == GLFW_RELEASE)
    {
        if (r_mouse_state_ == GLFW_PRESS)
        {
            RightMouseRelease();
        }
    }
    r_mouse_state_ = current_r_mouse_state;
}

void MainFrame::LeftMouseClick(float x, float y)
{
    if (modeling_state_ == OBJ_SUBDIVIDE)
    {
        glm::vec3 p_world = Screen2World(x, y);
        glm::vec3 cam_pos = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
        mesh_.SubdivideFace(cam_pos, glm::normalize(p_world - cam_pos));
    }
    else if (modeling_state_ == OBJ_EXTRUDE)
    {
        glm::vec3 p_world = Screen2World(x, y);
        glm::vec3 cam_pos = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
        mesh_.GenExtrudeFace(cam_pos, glm::normalize(p_world - cam_pos));
    }
}

void MainFrame::LeftMouseRelease()
{
    mesh_.CommitTransform();
}

void MainFrame::RightMouseClick(float x, float y)
{
    return;
}

void MainFrame::RightMouseMove(float start_x, float start_y, float curr_x, float curr_y)
{
    glm::vec2 s_start(start_x, start_y);
    glm::vec2 s_cur(curr_x, curr_y);
    glm::vec2 V = s_cur - s_start;
    glm::vec2 A = glm::vec2(-V.y, V.x);
    glm::vec3 rot_axis = glm::normalize(Screen2World(A + s_start) - Screen2World(s_start));
    glm::mat4x4 rot_mat = glm::rotate(glm::mat4x4(1.f), 0.007f * glm::length(A), rot_axis);
    camera_.ApplyTransform(rot_mat);
}

void MainFrame::RightMouseRelease()
{
    return;
}

glm::vec3 MainFrame::Camera2World(const glm::vec3 &x, float w)
{
    return glm::vec3(camera_.view_mat_inv_ * glm::vec4(x, w));
}

glm::vec3 MainFrame::World2Camera(const glm::vec3 &x, float w)
{
    return glm::vec3(camera_.view_mat_ * glm::vec4(x, w));
}

glm::vec3 MainFrame::Screen2World(const glm::vec2 &v, float depth)
{
    float x = v.x / SCR_WIDTH * 2.f - 1.f;
    float y = v.y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 1.f);
    v_camera = v_camera * depth;
    glm::vec4 v_world = camera_.view_mat_inv_ * v_camera;
    return glm::vec3(v_world);
}

glm::vec3 MainFrame::Screen2World(float scr_x, float scr_y, float camera_z)
{
    float x = scr_x / SCR_WIDTH * 2.f - 1.f;
    float y = scr_y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 1.f);
    v_camera = v_camera * -camera_z;
    glm::vec4 v_world = camera_.view_mat_inv_ * v_camera;
    return glm::vec3(v_world);
}

std::tuple<glm::vec3, glm::vec3> MainFrame::Screen2WorldRay(float scr_x, float scr_y)
{
    float x = scr_x / SCR_WIDTH * 2.f - 1.f;
    float y = scr_y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec3 o = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 0.f);
    glm::vec3 v = camera_.view_mat_inv_ * v_camera;
    return std::make_tuple(o, v);
}

void MainFrame::DrawScene()
{
    // Draw mesh
    mesh_.Draw();

    VisualizeWorldSpace();
}