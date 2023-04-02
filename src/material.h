#ifndef MATERIAL_H
#define MATERIAL_H

#include "rtweekend.h"

#include "texture.h"
#include "hittable.h"

enum Scatter_type
{
    REFRACT,
    REFLECT,
    DIFFUSE,
};

struct scatter_record
{
    ray scatter_ray;
    Scatter_type scatter_type;
    color attenuation;
    double brdf;
    double pdf;
};

// class material
// {
// public:
//     virtual color emitted(
//         const ray &r_in, const hit_record &rec, double u, double v, const point3 &p) const
//     {
//         return color(0, 0, 0);
//     }

//     virtual bool scatter(
//         const ray &r_in, const hit_record &rec, scatter_record &srec) const
//     {
//         return false;
//     }

//     virtual double scattering_pdf(
//         const ray &r_in, const hit_record &rec, const ray &scattered) const
//     {
//         return 0;
//     }
// };

class material
{
public:
    vec3 kd, ks, ke, tr;
    double ni, ns;
    double threshold;
    double preflec, pdiffus;
    shared_ptr<image_texture> map_kd;
    bool has_texture;
    double R0;
    material(const vec3 &_kd, const vec3 &_ks, const vec3 &_ke,
             const double &_ni, const double &_ns, const vec3 &_tr, bool has_texure_, shared_ptr<image_texture> _mat_texture)
        : kd(_kd), ks(_ks), ke(_ke), ni(_ni), ns(_ns), tr(_tr), has_texture(has_texure_), map_kd(_mat_texture)
    {
        if (ns <= 1)
        {
            pdiffus = 1.f;
        }
        else
        {
            pdiffus = 1.0f / (std::log10(ns) + 1.0f);
        }
        preflec = 1.0 - pdiffus;
        if (hasRefraction())
        {
            R0 = std::powf((1 - ni) / (1 + ni), 2);
            // prefrac = ni / (ni + 5.0);
        }
        else
        {
            R0 = 0;
        }
    }

    double random_scatter_type(Scatter_type &st, double cosi, double &pdf)
    {
        double Ri = R0 + (1 - R0) * std::powf(1 - cosi, 5);
        if (hasRefraction())
        {
            double p = random_double();
            if (p > Ri)
            {
                st = REFRACT;
                pdf = (1 - Ri);
            }
            else
            {
                double p1 = random_double();
                if (p1 < pdiffus)
                {
                    st = DIFFUSE;
                    pdf = pdiffus * Ri;
                }
                else
                {
                    st = REFLECT;
                    pdf = preflec * Ri;
                }
            }
            return Ri;
        }
        else
        {
            double p = random_double();
            if (p < pdiffus)
            {
                st = DIFFUSE;
                pdf = pdiffus;
            }
            else
            {
                st = REFLECT;
                pdf = preflec;
            }
            return 1.0;
        }
    }

    double getThreshold() { return threshold; }

    bool hasemission() const
    {
        return !ke.near_zero();
    }

    vec3 eval(vec3 wi, vec3 wo, vec3 n, Scatter_type s_type)
    {
        if (s_type == DIFFUSE)
        {
            return kd / PI;
        }
        else if (s_type == REFLECT)
        {
            vec3 refl = reflect(wo, n);
            double cosine = fmax(dot(wi, refl), 0.0f);
            return ks * std::pow(cosine, ns) * (+1.0f) / (PI * 2.0f);
        }
        else
        {
            return vec3();
        }
    }

    bool hasTexture()
    {
        return has_texture;
    }

    color getEmission() const
    {
        return ke;
    }

    color emitted(
        const ray &r_in, const hit_record &rec, double u, double v, const point3 &p) const
    {
        if (hasemission())
        {
            return ke;
        }
        else
        {
            return vec3(0, 0, 0);
        }
    }

    bool hasRefraction()
    {
        return ni > 1.0;
    }

    vec3 brdf(const ray &r, const hit_record &rec, const scatter_record &srec)
    {
        vec3 wi = -r.dir, wo = srec.scatter_ray.dir;
        vec3 h = normalize(wi + wo);
        double cos_nh = fmax(0.0, dot(rec.normal, h));
        vec3 diffuse = kd;
        vec3 specular = ks;
        if (has_texture)
        {
            diffuse = diffuse * map_kd->value(rec.uv.x(), rec.uv.y(), rec.p);
        }
        vec3 temp_color = diffuse + specular;
        // return diffuse / PI + specular * (ns + 1) * std::powf(dot(rec.normal, h), ns) / (8 * PI * std::powf(dot(wo, h), 3));
        specular *= std::pow(cos_nh, ns);
        return (diffuse / PI) + specular * ((ns + 2.0f) * (ns + 4.0f) / (8.0f * PI * (ns + std::pow(2.0, -ns / 2.0))));
        // +specular *((ns + 2.0f) * (ns + 4.0f) / (8.0f * M_PI * (ns + std::pow(2.0, -ns / 2.0))));
    }

    void scatter_(
        const ray &r_in, const hit_record &rec, scatter_record &srec)
    {
        Scatter_type st;
        double pdf;
        double cosi = dot(rec.normal, -r_in.dir);
        double Ri = random_scatter_type(st, cosi, pdf);
        // srec.pdf = pdf;
        srec.scatter_ray.orig = rec.p;
        srec.pdf = pdf;
        vec3 wi = -r_in.dir;
        vec3 wo;
        double inv_ni = 1 / ni;
        double p = pdiffus / (pdiffus + preflec);
        if (hasRefraction())
        {
            // std::cout << "refra";
            vec3 temp_n = rec.normal;
            if (dot(wi, rec.normal) < 0)
            {
                temp_n = -temp_n;
                inv_ni = ni;
            }
            double cos0 = dot(wi, temp_n);
            double cos1 = 1 - inv_ni * inv_ni * (1 - cos0 * cos0);
            if (cos1 > 0)
            {
                wo = temp_n * (inv_ni * cos0 - std::sqrt(cos1)) - wi * inv_ni;
            }
            else
            {
                wo = temp_n * 2 * ni - wi;
            }
            srec.scatter_ray.dir = wo;
            srec.pdf = 1.0;
            srec.scatter_type = REFRACT;
            // srec.scatter_ray.dir = random_refract(-r_in.dir, rec.normal, ni);
            // srec.pdf = prefrac;
            // srec.scatter_type = REFRACT;
        }
        // else if (st == DIFFUSE)
        else if (random_double() < pdiffus)
        {
            srec.scatter_ray.dir = random_hemisphere(-r_in.dir, rec.normal, ns, srec.pdf, p);
            // srec.pdf *= Ri;
            srec.scatter_type = DIFFUSE;
        }
        else
        {
            srec.scatter_ray.dir = random_hemisphere_specular(-r_in.dir, rec.normal, ns, srec.pdf, p);
            // srec.pdf *= Ri;
            srec.scatter_type = REFLECT;
        }
    }

    void scatter(
        const ray &r_in, const hit_record &rec, scatter_record &srec)
    {
        srec.scatter_ray.orig = rec.p;
        srec.pdf = 1.0;
        vec3 wi = -r_in.dir;
        double p;
        if (ks.near_zero())
        {
            p = 1;
        }
        else
        {
            if (ns <= 1)
            {
                p = 1;
            }
            else
            {
                p = 1.0 / (std::log10(ns) + 1.0);
            }
        }
        double refract_p = .0;
        bool is_refract = false;
        // if (hasRefraction())
        if (hasRefraction())
        {
            // std::cout << "fract" << std::endl;
            refract_p = 1.0;
            vec3 temp_n = rec.normal;
            double inv_ni = 1.0 / ni;
            if (dot(rec.normal, wi) < 0)
            {
                temp_n = -rec.normal;
                inv_ni = ni;
                refract_p = 1.0;
            }
            if (random_double() <= refract_p)
            {
                srec.pdf *= refract_p;
                double nl = dot(wi, temp_n);
                double temp = 1 - inv_ni * inv_ni * (1 - nl * nl);
                vec3 wo;
                if (temp > 0)
                {
                    wo = temp_n * (inv_ni * nl - sqrt(temp)) - wi * inv_ni;
                }
                else
                {
                    wo = temp_n * (2 * ni) - wi;
                }
                is_refract = true;
                srec.scatter_type = REFRACT;
                srec.scatter_ray.dir = wo;
                return;
            }
        }

        // diffuse
        // if (random_double() < p)
        // // if (1)
        // {
        //     // std::cout << "diffuse" << std::endl;
        //     double costheta;
        //     vec3 wo = random_hemisphere(rec.normal, 1.0, costheta);
        //     // srec.pdf = p * sintheta / TWO_PI;
        //     srec.pdf = p * costheta / PI;
        //     if (ns > 1)
        //     {
        //         // std::cout << "wdnmd";
        //         vec3 half = normalize(wi + wo);
        //         double temp_cos1 = dot(half, rec.normal);
        //         double temp_cos2 = dot(half, wi);
        //         srec.pdf += (1 - p) * (ns + 2) / (2 * PI) * std::pow(temp_cos1, ns + 1) / (4 * temp_cos2);
        //     }
        //     srec.pdf *= (1 - refract_p);
        //     srec.pdf = fmax(0, srec.pdf);
        //     // srec.pdf = 1.0;
        //     srec.scatter_ray.dir = wo;
        //     // srec.scatter_ray.dir = reflect(r_in.dir, rec.normal);
        //     srec.scatter_type = DIFFUSE;
        //     return;
        // }
        // else
        // {
        //     // std::cout << "specular" << std::endl;
        //     double costheta;
        //     vec3 local = random_hemisphere(rec.normal, ns, costheta);
        //     double temp_cos = dot(local, wi);
        //     vec3 wo = local * 2 * temp_cos - wi;

        //     srec.pdf = (1 - p) * (ns + 2) / (2 * PI) * std::pow(costheta, ns + 1);
        //     srec.pdf /= (4 * temp_cos);
        //     if (dot(rec.normal, wo) <= 0)
        //     {
        //         srec.pdf = 0;
        //     }
        //     else
        //     {
        //         srec.pdf += p * dot(rec.normal, wo) / PI;
        //     }
        //     srec.pdf *= (1 - refract_p);
        //     srec.scatter_ray.dir = wo;
        //     srec.scatter_type = REFLECT;
        //     return;
        // }
    }
};

// class lambertian : public material
// {
// public:
//     lambertian(const color &a) : albedo(make_shared<solid_color>(a)) {}
//     lambertian(shared_ptr<texture> a) : albedo(a) {}

//     virtual bool scatter(
//         const ray &r_in, const hit_record &rec, scatter_record &srec) const override
//     {
//         srec.is_specular = false;
//         srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
//         srec.pdf_ptr = make_shared<cosine_pdf>(rec.normal);
//         return true;
//     }

//     double scattering_pdf(
//         const ray &r_in, const hit_record &rec, const ray &scattered) const override
//     {
//         auto cosine = dot(rec.normal, unit_vector(scattered.direction()));
//         return cosine < 0 ? 0 : cosine / pi;
//     }

// public:
//     shared_ptr<texture> albedo;
// };

// class metal : public material
// {
// public:
//     metal(const color &a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

//     virtual bool scatter(
//         const ray &r_in, const hit_record &rec, scatter_record &srec) const override
//     {
//         vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
//         srec.specular_ray =
//             ray(rec.p, reflected + fuzz * random_in_unit_sphere(), r_in.time());
//         srec.attenuation = albedo;
//         srec.is_specular = true;
//         srec.pdf_ptr = nullptr;
//         return true;
//     }

// public:
//     color albedo;
//     double fuzz;
// };

// class dielectric : public material
// {
// public:
//     dielectric(double index_of_refraction) : ir(index_of_refraction) {}

//     virtual bool scatter(
//         const ray &r_in, const hit_record &rec, scatter_record &srec) const override
//     {
//         srec.is_specular = true;
//         srec.pdf_ptr = nullptr;
//         srec.attenuation = color(1.0, 1.0, 1.0);
//         double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

//         vec3 unit_direction = unit_vector(r_in.direction());
//         double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
//         double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

//         bool cannot_refract = refraction_ratio * sin_theta > 1.0;
//         vec3 direction;

//         if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
//             direction = reflect(unit_direction, rec.normal);
//         else
//             direction = refract(unit_direction, rec.normal, refraction_ratio);

//         srec.specular_ray = ray(rec.p, direction, r_in.time());
//         return true;
//     }

// public:
//     double ir; // Index of Refraction

// private:
//     static double reflectance(double cosine, double ref_idx)
//     {
//         // Use Schlick's approximation for reflectance.
//         auto r0 = (1 - ref_idx) / (1 + ref_idx);
//         r0 = r0 * r0;
//         return r0 + (1 - r0) * pow((1 - cosine), 5);
//     }
// };

// class diffuse_light : public material
// {
// public:
//     diffuse_light(shared_ptr<texture> a) : emit(a) {}
//     diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}

//     virtual color emitted(
//         const ray &r_in, const hit_record &rec, double u, double v, const point3 &p) const override
//     {
//         if (!rec.front_face)
//             return color(0, 0, 0);
//         return emit->value(u, v, p);
//     }

// public:
//     shared_ptr<texture> emit;
// };

// class isotropic : public material
// {
// public:
//     isotropic(color c) : albedo(make_shared<solid_color>(c)) {}
//     isotropic(shared_ptr<texture> a) : albedo(a) {}

// #if 0
//         // Issue #669
//         // This method doesn't match the signature in the base `material` class, so this one's
//         // never actually called. Disabling this definition until we sort this out.

//         virtual bool scatter(
//             const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
//         ) const override {
//             scattered = ray(rec.p, random_in_unit_sphere(), r_in.time());
//             attenuation = albedo->value(rec.u, rec.v, rec.p);
//             return true;
//         }
// #endif

// public:
//     shared_ptr<texture> albedo;
// };

#endif
