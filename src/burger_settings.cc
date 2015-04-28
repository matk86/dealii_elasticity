#include <burger.h>

namespace burger
{
  using namespace dealii;

  //LameLambda::LameLambda
  template <int dim>
  LameLambda<dim>::LameLambda () : Function<dim> ()
  {}

  //LameLambda::value
  template <int dim>
  inline
  double 
  LameLambda<dim>::value (const Point<dim> &p, 
			  const unsigned int ) const
  {
    //if (p[dim-1] > 0)
    if (p[dim-1] > 0.55 || p[dim-1] < -0.55)
      return 423.33;//ZnS
    //return 403.67;//CdS
    else
      return 145.67;//CuS
  }

  // LameLambda::value_list
  template <int dim>
  void 
  LameLambda<dim>::value_list (const std::vector<Point<dim> > &points,
			       std::vector<double>  &value_list, 
			       const unsigned int ) const
  {
    const unsigned int n_points = points.size();
    for (unsigned int p=0; p<n_points; ++p)
      value_list[p] = value(points[p]);
  }

  // LameMu::LameMu
  template <int dim>
  LameMu<dim>::LameMu () : Function<dim> ()
  {}

  //LameMu::value
  template <int dim>
  inline
  double 
  LameMu<dim>::value (const Point<dim> &p, const unsigned int) const
  {
    //if (p[dim-1] > 0)
    if (p[dim-1] > 0.55 || p[dim-1] < -0.55)
      return 290.0;//ZnS
    //return 166.33;//CdS
    else
      return 363.67;//CuS
  }

  // LameMu::value_list
  template <int dim>
  void 
  LameMu<dim>::value_list (const std::vector<Point<dim> > &points, 
			   std::vector<double>   &value_list, 
			   const unsigned int) const
  {
    const unsigned int n_points = points.size();
    for (unsigned int p=0; p<n_points; ++p)
      value_list[p] = value(points[p]);
  }

  //constructor
  template <int dim>
  Temperature<dim>::Temperature () : Function<dim> ()
  {}

  //Temperature<dim>::value
  template <int dim>
  inline
  double 
  Temperature<dim>::value (const Point<dim> &p, 
				  const unsigned int ) const
  {
    //if ((p[2] > -0.1 && p[2] < 0.1))
    //if (p[dim-1] > 0)
    if (p[dim-1] > 0.55 || p[dim-1] < -0.55)
      return -0.01; // strain in ZnS realtive to CuS
    else
      return 0;
  }

  //Temperature::value_list
  template <int dim>
  void 
  Temperature<dim>::value_list(const std::vector<Point<dim> > &points,
			       std::vector<double>  &value_list, 
			       const unsigned int ) const
  {
    const unsigned int n_points = points.size();
    for (unsigned int p=0; p<n_points; ++p)
      value_list[p] = value(points[p]);
  }

  template class LameLambda<3>;
  template class LameMu<3>;
  template class Temperature<3>;

}
