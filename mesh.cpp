#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_ply_points.h>
#include <utility>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Timer.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/jet_smooth_point_set.h>
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<unsigned char, 3> Color;
// Point with normal, color and intensity
typedef CGAL::cpp11::tuple<Point, Vector, Color, int> PNCI;
typedef CGAL::Nth_of_tuple_property_map<0, PNCI> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNCI> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNCI> Color_map;
typedef CGAL::Nth_of_tuple_property_map<3, PNCI> Intensity_map;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Scale_space_surface_reconstruction_3<Kernel> Reconstruction;


namespace CGAL {
    template< class F >
    struct Output_rep< ::Color, F > {
        const ::Color& c;
        static const bool is_specialized = true;
        Output_rep (const ::Color& c) : c(c)
        { }
        std::ostream& operator() (std::ostream& out) const
        {
            if (is_ascii(out))
                out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]) << " " << int(c[3]);
            else
                out.write(reinterpret_cast<const char*>(&c), sizeof(c));
            return out;
        }
    };
}


int main(int argc, char*argv[])
{
    
    const char* fname = (argc>1) ? argv[1] : "pills.ply";
    // Reads a .ply point set file with normal vectors and colors
    std::vector<PNCI> points; // store points
   // store points
    std::ifstream in(fname);
    if (!in ||
        !CGAL::read_ply_points_with_properties
        (in,
         std::back_inserter (points), // output
         CGAL::make_ply_point_reader (Point_map())// generates a property handler to read 3D points
         ))
    {
        std::cerr << "Error: cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }
    
    //copy the input point cloud vector into a different format
    std::vector<Point> points2(points.size());

            for (std::size_t i = 0; i < points.size(); ++ i)
            {
                const Point& p = get<0>(points[i]);
                Point &p2 = points2[i];
                Point p1(p.x(),p.y(),p.z());
                p2 = p1;
            }
    //thin the point cloud
    CGAL::Timer task_timer; task_timer.start();
    // simplification by clustering using erase-remove idiom
    points2.erase (CGAL::hierarchy_simplify_point_set (points2,
                                                      CGAL::parameters::size(100). // Max cluster size
                                                      maximum_variation(0.01)), // Max surface variation
                  points2.end ());

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::vector<PNCI> points3(points2.size());

    // write to output vector
    for (std::size_t i = 0; i < points3.size(); ++ i)
    {
        Point& p = get<0>(points3[i]);
        Point &p2 = points2[i];
        Point p1(p2.x(),p2.y(),p2.z());
        p = p1;
    }

    
    std::cout << points3.size () << " point(s) kept, computed in "
    << task_timer.time() << " seconds, "
    << (memory>>20) << " Mib allocated." << std::endl;
    
    std::ofstream f("thinnedPLY.ply");
    CGAL::set_ascii_mode(f); // The PLY file will be written in the binary format

        if(!CGAL::write_ply_points_with_properties
        (f, points3,
         CGAL::make_ply_point_writer (Point_map()))){
            std::cerr<<"errrrr";
        }

        std::cerr << "done: " << points2.size() << " points." << std::endl;
        std::cerr << "Reconstruction ";
        CGAL::Timer t;
        t.start();
        // Construct the mesh in a scale space.
        Reconstruction reconstruct (points2.begin(), points2.end());
        reconstruct.increase_scale(4);
        reconstruct.reconstruct_surface();
        std::cerr << "done in " << t.time() << " sec." << std::endl;
        t.reset();
        std::ofstream out ("outMesh.off");
        out << reconstruct;
        std::cerr << "Writing result in " << t.time() << " sec." << std::endl;
        std::cerr << "Done." << std::endl;
        return EXIT_SUCCESS;
}




