// Raphaele Michelle Guillemot
// COMP3271 Programming Assignment 3
// Hittable.cpp

#include "Hittable.h"

// Sphere
bool Sphere::Hit(const Ray &ray, HitRecord *hit_record) const
{
    bool ret = false;

    // get relative center position
    glm::vec3 p_ = ray.o - o_;
    float a = glm::dot(ray.d, ray.d);
    float b = 2 * glm::dot(ray.d, p_);
    float c = glm::dot(p_, p_) - (r_ * r_);
    float discriminant = (b * b) - (4 * a * c);
    float t = 0.0f;

    if (discriminant == 0.0f)
    {
        t = -b / (2 * a);
        if (t > 0)
            ret = true;
    }
    else if (discriminant > 1e-5f)
    {
        float t1 = (-b + sqrt(discriminant)) / (2 * a);
        float t2 = (-b - sqrt(discriminant)) / (2 * a);
        if (t1 < t2)
            t = t1;
        else
            t = t2;
        if (t > 0)
            ret = true;
    }

    if (ret)
    {
        hit_record->position = ray.o + t * ray.d;
        hit_record->normal = glm::normalize(hit_record->position - o_);
        hit_record->distance = t;
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal);
        hit_record->material = material_;
    }

    return ret;
}

// Quadric
bool Quadric::Hit(const Ray &ray, HitRecord *hit_record) const
{
    bool ret = false;
    glm::vec4 O(ray.o, 1.0f);
    glm::vec4 D(ray.d, 0.0f);
    float a = glm::dot(D, A_ * D);
    float b = glm::dot(O, 2.0f * A_ * D);
    float c = glm::dot(O, A_ * O);
    float discriminant = (b * b) - (4 * a * c);
    float t = 0.0f;

    if (discriminant == 0.0f)
    {
        t = -b / (2 * a);
        if (t > 1e-5)
            ret = true;
    }
    else if (discriminant > 1e-5f)
    {
        float t1 = (-b + sqrt(discriminant)) / (2 * a);
        float t2 = (-b - sqrt(discriminant)) / (2 * a);
        if (t1 < t2)
            t = t1;
        else
            t = t2;
        if (t > 1e-5f)
            ret = true;
    }

    if (ret)
    {
        hit_record->position = ray.o + t * ray.d;
        hit_record->normal = glm::normalize((A_ + glm::transpose(A_)) * glm::vec4(hit_record->position, 1.0f));
        hit_record->distance = t;
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal);
        hit_record->material = material_;
    }
    return ret;
}

// Triangle
bool Triangle::Hit(const Ray &ray, HitRecord *hit_record) const
{
    bool ret = false;
    glm::vec3 n = glm::cross(b_ - a_, c_ - a_);
    if (glm::dot(ray.d, n) != 0.0f)
    {
        // d is the scalar
        float d = glm::dot(n, a_);
        // t derived from ray equation
        float t = (d - glm::dot(n, ray.o)) / glm::dot(n, ray.d);
        // float t = glm::dot(a_ - ray.o, n) / glm::dot(ray.d, n);
        if (t > 1e-5f)
        {
            glm::vec3 p = ray.o + t * ray.d;
            glm::vec3 oa = p - a_;
            glm::vec3 ob = p - b_;
            glm::vec3 oc = p - c_;
            glm::vec3 cross_oab = glm::cross(oa, ob);
            glm::vec3 cross_obc = glm::cross(ob, oc);
            glm::vec3 cross_oca = glm::cross(oc, oa);
            // all same direction -> in triangle
            if (glm::dot(cross_oab, cross_obc) > 0 && glm::dot(cross_obc, cross_oca) > 0 && glm::dot(cross_oca, cross_oab) > 0)
            {
                // calculate normal
                if (phong_interpolation_)
                {
                    float oab = glm::length(glm::cross(oa, ob)) / 2;
                    float obc = glm::length(glm::cross(ob, oc)) / 2;
                    float oca = glm::length(glm::cross(oc, oa)) / 2;
                    float area = glm::length(n) / 2;

                    float alphaA = obc / area;
                    float alphaB = oca / area;
                    float alphaC = oab / area;

                    hit_record->normal = alphaA * n_a_ + alphaB * n_b_ + alphaC * n_c_;
                }
                else
                    hit_record->normal = glm::normalize(n);

                ret = true;
                hit_record->position = p;
                hit_record->reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal);
                hit_record->distance = t;
                hit_record->in_direction = ray.d;
            }
        }
    }
    return ret;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray &ray, HitRecord *hit_record) const
{
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret)
    {
        hit_record->material = material_;
    }
    return ret;
}

// Mesh
Mesh::Mesh(const std::string &file_path,
           const Material &material,
           bool phong_interpolation) : ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation)
{
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++)
    {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto &face : f_ind_)
    {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto &vertex_normal : vertex_normals_)
    {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto &face : f_ind_)
    {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
                                vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
                                phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min(1e5f, 1e5f, 1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto &vertex : vertices_)
    {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++)
    {
        InsertFace(root_, i);
    }
}

bool Mesh::Hit(const Ray &ray, HitRecord *hit_record) const
{
    const bool brute_force = false;
    if (brute_force)
    {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto &triangle : triangles_)
        {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record))
            {
                if (curr_hit_record.distance < min_dist)
                {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f)
        {
            hit_record->material = material_;
            return true;
        }
        return false;
    }
    else
    {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret)
        {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t> &face, const Point &bbox_min, const Point &bbox_max) const
{
    for (size_t idx : face)
    {
        const auto &pt = vertices_[idx];
        for (int i = 0; i < 3; i++)
        {
            if (pt[i] < bbox_min[i] + 1e-6f)
                return false;
            if (pt[i] > bbox_max[i] - 1e-6f)
                return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray &ray, const Point &bbox_min, const Point &bbox_max) const
{
    float t_min = -1e5f;
    float t_max = 1e5f;

    for (int i = 0; i < 3; i++)
    {
        if (glm::abs(ray.d[i]) < 1e-6f)
        {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f)
            {
                t_min = 1e5f;
                t_max = -1e5f;
            }
        }
        else
        {
            if (ray.d[i] > 0.f)
            {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else
            {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode *u, size_t face_idx)
{
    const Point &bbox_min = u->bbox_min;
    const Point &bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++)
    {
        for (size_t b = 0; b < 2; b++)
        {
            for (size_t c = 0; c < 2; c++)
            {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max))
                {
                    if (u->childs[child_idx] == nullptr)
                    {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode *child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs)
    {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode *u, const Ray &ray, HitRecord *hit_record) const
{
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max))
    {
        return false;
    }
    float distance = 1e5f;
    for (const auto &face_idx : u->face_index)
    {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < distance)
            {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto &child : u->childs)
    {
        if (child == nullptr)
        {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < distance)
            {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}

// Hittable list
void HittableList::PushHittable(const Hittable &hittable)
{
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray &ray, HitRecord *hit_record) const
{
    float min_dist = 1e5f;
    for (const auto &hittable : hittable_list_)
    {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < min_dist)
            {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}