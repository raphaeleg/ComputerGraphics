// Raphaele Michelle Guillemot
// COMP3271 Programming Assignment 3
// main.cpp

#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray &ray,
               const std::vector<LightSource> &light_sources,
               const Hittable &scene,
               int trace_depth);


Color Shade(const std::vector<LightSource> &light_sources,
            const Hittable &hittable_collection,
            const HitRecord &hit_record,
            int trace_depth)
{
    Color color(0.f, 0.f, 0.f);

    // ************
    // Ambient term
    // ************
    color += hit_record.material.k_a * hit_record.material.ambient;

    // ************
    // Diffuse and Specular
    // ************

    // for each light source
    for (int i = 0; i < light_sources.size(); i++)
    {
        // shadow ray = ray from intersect point to light source
        Ray shadowRay = Ray(hit_record.position, light_sources[i].position - hit_record.position);
        glm::vec3 shadowRayN = glm::normalize(shadowRay.d);
        // if N.dot(Li) > 0
        if (glm::dot(hit_record.normal, shadowRayN) > 0)
        {
            // if no intersection between shadow ray and scene
            HitRecord shadowHitRecord;
            if (!hittable_collection.Hit(shadowRay, &shadowHitRecord))
            {
                Color diffuseTerm = hit_record.material.k_d * hit_record.material.diffuse * glm::dot(hit_record.normal, shadowRayN);
                
                // R = direction of the reflection of shadow ray
                // normal * (2 * (normal * incident)) - incident
                // glm::vec3 R = shadowRay.d - (hit_record.normal) * (2 * glm::dot(hit_record.normal, shadowRay.d));
                glm::vec3 R;
                R[0] = 2 * hit_record.normal[0] * (hit_record.normal[0] * shadowRay.d[0] + hit_record.normal[1] * shadowRay.d[1] + hit_record.normal[2] * shadowRay.d[2]) - shadowRay.d[0];
                R[1] = 2 * hit_record.normal[1] * (hit_record.normal[0] * shadowRay.d[0] + hit_record.normal[1] * shadowRay.d[1] + hit_record.normal[2] * shadowRay.d[2]) - shadowRay.d[1];
                R[2] = 2 * hit_record.normal[2] * (hit_record.normal[0] * shadowRay.d[0] + hit_record.normal[1] * shadowRay.d[1] + hit_record.normal[2] * shadowRay.d[2]) - shadowRay.d[2];
                R = glm::normalize(R);
                
                // V = reversed direction of the shoot in ray
                glm::vec3 V = -(hit_record.in_direction);
                // V = glm::normalize(V);

                // clamp
                float specularConst = glm::dot(R, V);
                if (specularConst <= 1e-5f)
                    specularConst = 0.0f;

                Color specularTerm = hit_record.material.k_s * hit_record.material.specular * pow(specularConst, hit_record.material.sh);
                color += light_sources[i].intensity * (diffuseTerm + specularTerm);
            }
        }
    }

    // ************
    // Relected Rays
    // ************
    if (trace_depth < kMaxTraceDepth)
    {
        if (hit_record.material.k_s > 1e-5f)
        {                                                                                                // i.e., k_s > 0
            Ray reflectedRay = Ray(hit_record.position, hit_record.reflection);                          // ray in reflection direction from intersection;
            Color r_color = TraceRay(reflectedRay, light_sources, hittable_collection, trace_depth + 1); 
            color += hit_record.material.k_s * r_color; // scale r_color by reflectance and add to color;
        }
    }

    // Clamp color to 1
    for (int i = 0; i < 3; i++)
    {
        if (color[i] > 1.0f)
            color[i] = 1.0f;
    }

    return color;
}

Color TraceRay(const Ray &ray,
               const std::vector<LightSource> &light_sources,
               const Hittable &hittable_collection,
               int trace_depth)
{
    // IF hittable_collection.Hit(...) Then Call shade function to calculate the color of the intersection point. 􏰕Return this color.
    // ELSE Return the background color(0,0,0).
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);

    if (hittable_collection.Hit(ray, &record))
        color = Shade(light_sources, hittable_collection, record, trace_depth);

    return color;
}

int main()
{
    // TODO: Set your workdir (absolute path) here.
    const std::string work_dir("/Users/crushedsummers/Desktop/Graphics_PA3_Release/");

    // Construct scene
    Scene scene(work_dir, "scene/teapot.toml");
    const Camera &camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f)
            {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f)
                {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++)
            {
                image[idx + i] = (uint8_t)(glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f)
            {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "/outputs/myteapot.png");
}
