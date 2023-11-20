#pragma once
#include "Common.h"
#include "Ray.h"


// Simple material
// Reference: https://en.wikipedia.org/wiki/Phong_reflection_model
struct Material {
    Material() = default;
    Material(const Color& ambient,
             const Color& diffuse,
             const Color& specular,
             float k_a,
             float k_d,
             float k_s,
             float sh) :
             ambient(ambient), diffuse(diffuse), specular(specular), k_a(k_a), k_d(k_d), k_s(k_s), sh(sh) {
    }

    Color ambient, diffuse, specular;
    float k_a, k_d, k_s;
    float sh = 1.f;             // shiness
};
