#include <cmath>

#include <deal.II/base/tensor.h>

#include <ale.h>
#include <utility.h>

namespace ale
{

using namespace dealii;

//LameLambda::LameLambda
template <int dim>
LameLambda<dim>::LameLambda (double val) : Function<dim> ()
{
	ALE_LAMBDA = val;
}

// LameMu::LameMu
template <int dim>
LameMu<dim>::LameMu (double val) : Function<dim> ()
{
	ALE_MU = val;
}

// BodyForce::BodyForce
template <int dim>
BodyForce<dim>::BodyForce(Point<dim> p1, Point<dim> p2, Point<dim> p3,
		std::vector<double> hkl,
			  double line_tol, double fmag, double fratio)
		{
	ALE_FORCE = fmag;
	ALE_FORCE_RATIO = fratio;
	ALE_TOL = line_tol;

	ALE_P1 = p1;   ALE_P2 = p2;   ALE_P3 = p3;

	// vector along line 1
	ALE_T1[0] =  ALE_P3[0] -  ALE_P1[0];
	ALE_T1[1] =  ALE_P3[1] -  ALE_P1[1];
	ALE_T1[2] =  ALE_P3[2] -  ALE_P1[2];

	// vector along line 2
	ALE_T2[0] =  ALE_P1[0] -  ALE_P2[0];
	ALE_T2[1] =  ALE_P1[1] -  ALE_P2[1];
	ALE_T2[2] =  ALE_P1[2] -  ALE_P2[2];

	// vector along line 3
	ALE_T3[0] =  ALE_P2[0] -  ALE_P3[0];
	ALE_T3[1] =  ALE_P2[1] -  ALE_P3[1];
	ALE_T3[2] =  ALE_P2[2] -  ALE_P3[2];

	// normal to the plane
	ALE_HKL[0] = hkl[0];  ALE_HKL[1] = hkl[1];  ALE_HKL[2] = hkl[2];
		}

//LameLambda::value
template <int dim>
inline
double
LameLambda<dim>::value (const Point<dim> &p,
		const unsigned int ) const
		{
	return ALE_LAMBDA;
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

//LameMu::value
template <int dim>
inline
double
LameMu<dim>::value (const Point<dim> &p, const unsigned int) const
{
	return ALE_MU;
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

// check whether the point is on the line defined by pbegin and pend
//points
template <int dim>
bool
BodyForce<dim>::isonline(const Point<dim> & p,
		const Point<dim> & pbegin,
		const Point<dim> & pend,
		double tol)
		{
	double c,x,y,z;
	// take care of lines parallel to xy, xz and yz planes first
	// parallel to yz
	if (abs(pend[0] - pbegin[0])<=1e-6)
	{
		c = (p[1] - pbegin[1]) / (pend[1] - pbegin[1]);
		x = pbegin[0];
		y = pbegin[1] + c * (pend[1] - pbegin[1]);
		z = pbegin[2] + c * (pend[2] - pbegin[2]);
	}
	// parallel to xz
	else if (abs(pend[1] - pbegin[1])<=1e-6)
	{
		c = (p[0] - pbegin[0]) / (pend[0] - pbegin[0]);
		x = pbegin[0] + c * (pend[0] - pbegin[0]);
		y = pbegin[1];
		z = pbegin[2] + c * (pend[2] - pbegin[2]);
	}
	// parallel to xy
	else if (abs(pend[2] - pbegin[2])<=1e-6)
	{
		c = (p[0] - pbegin[0]) / (pend[0] - pbegin[0]);
		x = pbegin[0] + c * (pend[0] - pbegin[0]);
		y = pbegin[1] + c * (pend[1] - pbegin[1]);
		z = pbegin[2] ;
	}
	else
	{
		c = (p[0] - pbegin[0]) / (pend[0] - pbegin[0]);
		x = p[0];
		y = pbegin[1] + c * (pend[1] - pbegin[1]);
		z = pbegin[2] + c * (pend[2] - pbegin[2]);
	}
	/*
	Tensor<1,3, double> t1, t2, t3;
	t1[0] = p[0] - pbegin[0];
	t1[1] = p[1] - pbegin[1];
	t1[2] = p[2] - pbegin[2];
	t2[0] = pend[0] - pbegin[0];
	t2[1] = pend[1] - pbegin[1];
	t2[2] = pend[2] - pbegin[2];
	cross_product<3, double>(t3,t1, t2);
	// sine of angle between vectors t1 and t2
	//double sint = t3.norm()/t1.norm()/t2.norm();
	 */

	if ( (abs(x-p[0]) <= tol) &&
			(abs(y-p[1]) <= tol) &&
			(abs(z-p[2]) <= tol) )
	{
		return true;
	}
	else
		return false;
		}


//returns the force vector at point p and normal to the line that
//contains p. checks only for 3 lines defined by p1, p2 and p3
// net force  =  direction normal to the line and the surface normal to 111
// and pointing outwards from the the 111 plane * 
//  surface_tension_111/surface_tension_100
// +  direction normal to the line and the surface normal to 001or 010 or 100
// and pointing outwards from the corresponding plane 
template <int dim>
std::vector<double>
BodyForce<dim>::get_force(const Point<dim> & p)
{
  Tensor<1,3, double> ftmp, ftmp1, ftmp2 ;
  Tensor<1,3, double> hkl001, hkl010, hkl100;
  hkl001[0] = 0.0;  hkl001[1] = 0.0;  hkl001[2] = 1.0;
  hkl010[0] = 0.0;  hkl010[1] = 1.0;  hkl010[2] = 0.0;
  hkl100[0] = 1.0;  hkl100[1] = 0.0;  hkl100[2] = 0.0;
  
  if (isonline(p,  ALE_P1,  ALE_P2,  ALE_TOL))
    {
      cross_product<3, double>(ftmp,  ALE_HKL,  ALE_T2);
      cross_product<3, double>(ftmp1,  ALE_T2, hkl010);
    }
  else if (isonline(p,  ALE_P2,  ALE_P3,  ALE_TOL))
    {
      cross_product<3, double>(ftmp,  ALE_HKL,  ALE_T3);
      cross_product<3, double>(ftmp1,  ALE_T3, hkl001);
    }
  else if (isonline(p,  ALE_P3,  ALE_P1,  ALE_TOL))
    {
      cross_product<3, double>(ftmp,  ALE_HKL,  ALE_T1);
      cross_product<3, double>(ftmp1,  ALE_T1, hkl100);
    }

    double ftmp_norm = sqrt(pow(ftmp[0],2) + pow(ftmp[1],2) + pow(ftmp[2],2));
    double ftmp1_norm = sqrt(pow(ftmp1[0],2) + pow(ftmp1[1],2) + pow(ftmp1[2],2));
  //  std::cout<< "normal : "<< ftmp[0] << " "
  //	       << ftmp[1] << " "
  //       << ftmp[2] << " " << std::endl;
  //std::cout << "the norm " << ftmp_norm << std::endl;
  // normalize
  if (abs(ftmp_norm) > 1e-6 && abs(ftmp1_norm) > 1e-6)
    {
      ftmp[0] = ftmp[0]/ftmp_norm;
      ftmp[1] = ftmp[1]/ftmp_norm;
      ftmp[2] = ftmp[2]/ftmp_norm;

      ftmp1[0] = ftmp1[0]/ftmp1_norm;
      ftmp1[1] = ftmp1[1]/ftmp1_norm;
      ftmp1[2] = ftmp1[2]/ftmp1_norm;

    }

  ftmp1 += ftmp * ALE_FORCE_RATIO;

  // force = normal * force magnitude
  // if the normal =0 then force = 0
  std::vector<double> force(3);
  force[0] = ftmp1[0] * ALE_FORCE;
  force[1] = ftmp1[1] * ALE_FORCE;
  force[2] = ftmp1[2] * ALE_FORCE;
  
  return force;
  
}


// list of force vectors for the given list of points
template<int dim>
void
BodyForce<dim>::force_value_list(const std::vector<Point<dim> > &points,
		std::vector<std::vector<double> >  &force_values)
		{
	const unsigned int n_points = points.size();
	for (unsigned int p=0; p<n_points; ++p)
	{
		force_values.push_back(get_force(points[p]));
	}

		}

// template instantiations for 3d
template class LameLambda<3>;
template class LameMu<3>;
template class BodyForce<3>;

}
