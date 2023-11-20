#pragma once
#include <string>
#include "Camera.h"
#include "Hittable.h"

struct LightSource {
    LightSource(const Point& position, float intensity) : position(position), intensity(intensity) {}
    Point position;
    float intensity;
};

// Scene = hittable objects + light sources + camera
class Scene {
public:
    // Constructor
    // Params:
    // file_path: path of the scene file
    Scene(const std::string& work_dir, const std::string& scane_file_name);

public:
    Camera camera_;
    std::vector<LightSource> light_sources_;
    HittableList hittable_collection_;

private:
    std::vector<std::unique_ptr<Hittable>> hittable_pool_;
    std::vector<std::unique_ptr<Material>> material_pool_;
};
