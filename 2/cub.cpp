#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t3");

  double lc = 1e-2;
  gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
  gmsh::model::geo::addPoint(.1, 0, 0, lc, 2);
  gmsh::model::geo::addPoint(0, .1, 0, lc, 3);
  gmsh::model::geo::addPoint(0, 0, .1, lc, 4);
  gmsh::model::geo::addPoint(0, .1, .1, lc, 5);
  gmsh::model::geo::addPoint(.1, 0, .1, lc, 6);
  gmsh::model::geo::addPoint(.1, .1, 0, lc, 7);
  gmsh::model::geo::addPoint(.1, .1, .1, lc, 8);

  gmsh::model::geo::addLine(1, 2, 1);
  gmsh::model::geo::addLine(1, 3, 2);
  gmsh::model::geo::addLine(1, 4, 3);
  gmsh::model::geo::addLine(2, 6, 4);
  gmsh::model::geo::addLine(2, 7, 5);
  gmsh::model::geo::addLine(4, 6, 6);
  gmsh::model::geo::addLine(4, 5, 7);
  gmsh::model::geo::addLine(3, 7, 8);
  gmsh::model::geo::addLine(3, 5, 9);
  gmsh::model::geo::addLine(8, 6, 10);
  gmsh::model::geo::addLine(8, 5, 11);
  gmsh::model::geo::addLine(8, 7, 12);

  gmsh::model::geo::addCurveLoop({1, 5, -8, -2}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);

  gmsh::model::geo::addCurveLoop({1, 4, -6, -3}, 2);
  gmsh::model::geo::addPlaneSurface({2}, 2);

  gmsh::model::geo::addCurveLoop({11, -7, 6, -10}, 3);
  gmsh::model::geo::addPlaneSurface({3}, 3);

  gmsh::model::geo::addCurveLoop({11, -9, 8, -12}, 4);
  gmsh::model::geo::addPlaneSurface({4}, 4);

  gmsh::model::geo::addCurveLoop({4, -10, 12, -5}, 5);
  gmsh::model::geo::addPlaneSurface({5}, 5);

  gmsh::model::geo::addCurveLoop({3, 7, -9, -2}, 6);
  gmsh::model::geo::addPlaneSurface({6}, 6);

  gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4, 5, 6}, 1);
  gmsh::model::geo::addVolume({1});

  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(3);

  gmsh::write("t3.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}

