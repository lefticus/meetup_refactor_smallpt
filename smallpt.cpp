#include <array>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2 

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm 
  double x{0}, y{0}, z{0};                  // position, also color (r,g,b) 
  constexpr Vec() noexcept = default;
  constexpr Vec(double x_, double y_=0, double z_=0) noexcept
    : x(x_), y(y_), z(z_)
  { 
  } 


  constexpr Vec operator+(const Vec b) const noexcept { 
    return Vec(x+b.x,y+b.y,z+b.z); 
  }

  constexpr Vec operator-(const Vec b) const noexcept { 
    return Vec(x-b.x,y-b.y,z-b.z); 
  }

  constexpr Vec operator*(double b) const noexcept { 
    return Vec(x*b,y*b,z*b); 
  }

  constexpr Vec mult(const Vec b) const noexcept { 
    return Vec(x*b.x,y*b.y,z*b.z); 
  } 

  constexpr Vec& norm() noexcept { 
    return *this = *this * (1/sqrt(x*x+y*y+z*z)); 
  }

  constexpr double dot(const Vec b) const noexcept { 
    return x*b.x+y*b.y+z*b.z; 
  }  

  // cross:
  constexpr Vec operator%(const Vec b) const noexcept {
    return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
  } 
}; 

struct Ray { 
  Vec o, d;          
  constexpr Ray(Vec o_, Vec d_) noexcept : o(o_), d(d_) {} 
};

enum class Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance() 

struct Sphere { 
  double rad;       // radius 
  Vec p, e, c;      // position, emission, color 
  Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 
  
  constexpr Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) noexcept
    : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) 
    {} 

  constexpr double intersect(const Ray &r) const noexcept { // returns distance, 0 if nohit 
    const Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
    constexpr const double eps=1e-4;
    const auto b=op.dot(r.d);
    auto det=b*b-op.dot(op)+rad*rad; 

    if (det<0) {
      return 0; 
    } else  {
      det=sqrt(det); 
    }

    if (const auto t = b-det; t > eps) {
      return t;
    } else if (const auto t = b+det; t > eps) {
      return t;
    } else {
      return 0;
    }
  } 
}; 

constexpr const Sphere spheres[] = {//Scene: radius, position, emission, color, material 
  Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),Refl_t::DIFF),//Left 
  Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),Refl_t::DIFF),//Rght 
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),Refl_t::DIFF),//Back 
  Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           Refl_t::DIFF),//Frnt 
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),Refl_t::DIFF),//Botm 
  Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),Refl_t::DIFF),//Top 
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, Refl_t::SPEC),//Mirr 
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, Refl_t::REFR),//Glas 
  Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), Refl_t::DIFF) //Lite 
}; 

constexpr double clamp(double x) noexcept
{ 
  return x<0 ? 0 : (x>1 ? 1 : x); 
} 

constexpr int toInt(double x) noexcept { 
  return int(pow(clamp(x),1/2.2)*255+.5); 
} 

constexpr bool intersect(const Ray &r, double &t, int &id) noexcept { 
  constexpr const auto n=sizeof(spheres)/sizeof(Sphere);
  constexpr const double inf=1e20;
  t=inf; 
  
  const auto close_enough = [](const auto a, const auto b)
  { 
    static_assert(std::is_same_v<std::common_type_t<decltype(a), decltype(b)>, double>);
    return std::abs(a-b) < std::numeric_limits<double>::epsilon();
  };

  for(auto i = n; i-- != 0; )  {
    
    if(const double d=spheres[i].intersect(r); !close_enough(d, 0.0) && d<t ){
      t=d;
      id=i;
    } 
  }
    
  return t < inf; 
} 

constexpr Vec radiance(const Ray &r, int depth, std::array<unsigned short, 3> &Xi)
{ 
  double t{};                               // distance to intersection 
  int id(0);                               // id of intersected object 
  
  if (!intersect(r, t, id)) return Vec(); // if miss, return black 
  
  const Sphere &obj = spheres[id];        // the hit object 
  
  const Vec x=r.o+r.d*t;
  const Vec n=(x-obj.p).norm();
  const Vec nl=n.dot(r.d)<0?n:n*-1;
  Vec f=obj.c; 
  const double p = (f.x>f.y && f.x>f.z) ? (f.x) : (f.y>f.z ? f.y : f.z); // max refl 
  
  if (++depth>5) {
    if (erand48(Xi.data())<p) {
      f=f*(1/p); 
    } else {
      return obj.e; //R.R. 
    }
  }
  
  if (obj.refl == Refl_t::DIFF){                  // Ideal DIFFUSE reflection 
    double r1=2*M_PI*erand48(Xi.data());
    double r2=erand48(Xi.data());
    double r2s=sqrt(r2); 
    const auto w=nl;
    const auto u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm();
    const auto v=w%u; 
    const Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
    return obj.e + f.mult(radiance(Ray(x,d),depth,Xi)); 
  } else if (obj.refl == Refl_t::SPEC) {
    // Ideal SPECULAR reflection 
    return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); 
  }

  const Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION 
  
  const bool into = n.dot(nl)>0;                // Ray from outside going in? 
  const double nc=1;
  const double nt=1.5;
  const double nnt= [&]{
    if (into) {
      return nc/nt;
    } else { 
      return nt/nc;
    }
  }();

  const double ddn=r.d.dot(nl);
  
  double cos2t{1-nnt*nnt*(1-ddn*ddn)}; 
  
  if (cos2t<0.0)    // Total internal reflection 
    return obj.e + f.mult(radiance(reflRay,depth,Xi)); 
  
  const Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
  
  const double a=nt-nc;
  const double b=nt+nc;
  const double R0=a*a/(b*b);
  const double c = 1-(into?-ddn:tdir.dot(n)); 
  const double Re=R0+(1-R0)*c*c*c*c*c;
  const double Tr=1-Re;
  const double P=.25+.5*Re;
  const double RP=Re/P;
  const double TP=Tr/(1-P); 
  
  return obj.e + f.mult(depth>2 ? (erand48(Xi.data())<P ?   // Russian roulette 
                                   radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) : 
                        radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 
} 

int main(int argc, char *argv[]){   
  const int w=160, h=120, samps = argc==2 ? atoi(argv[1])/4 : 60; // # samples 
  

  Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 
  Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h]; 

  for (int y=0; y<h; y++) {                      // Loop over image rows 
    
    fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); 
    
    std::array<unsigned short, 3> Xi{0,0,static_cast<unsigned short>(y*y*y)};
    for (int x=0; x<w; x++){   // Loop cols 
      const int i=(h-y-1)*w+x;
      for (int sy=0; sy<2; sy++){     // 2x2 subpixel rows 
        for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols 
          for (int s=0; s<samps; s++){ 
            const double r1=2*erand48(Xi.data());
            const double dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1); 
            const double r2=2*erand48(Xi.data());
            const double dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2); 
            Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) + 
              cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d; 
            r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps); 
          } // Camera rays are pushed ^^^^^ forward to start in interior 
          c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25; 
        } 
      }
    }
  } 

  FILE *f = fopen("/home/jason/image.ppm", "w");         // Write image to PPM file. 
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
  
  for (int i=0; i<w*h; ++i) { 
    fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
  }
   
  fclose(f);
} 
