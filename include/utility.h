#ifndef __utility_h
#define __utility_h

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>

#include <cmath>


using namespace dealii;

/*------------------------------------------------------------------
  Some utility functions
  -----------------------------------------------------------------*/

// return the elastic tensor corresponsind to the lame coefficents
template <int dim>
SymmetricTensor<4,dim>
  get_stress_strain_tensor (const double lambda, const double mu)
{
  SymmetricTensor<4,dim> tmp;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  tmp[i][j][k][l] = (((i==k) && (j==l) ? mu : 0.0) + 
			     ((i==l) && (j==k) ? mu : 0.0) + 
			     ((i==j) && (k==l) ? lambda : 0.0));
  return tmp;
}


// return the strain tensor
template <int dim>
inline
SymmetricTensor<2,dim>
  get_strain (const FEValues<dim> &fe_values, 
	      const unsigned int shape_func, 
	      const unsigned int q_point)
{
  SymmetricTensor<2,dim> tmp;
  for (unsigned int i=0; i<dim; ++i)
    tmp[i][i] = fe_values.shape_grad_component (shape_func,q_point,i)[i];
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i+1; j<dim; ++j)
      tmp[i][j] = (fe_values.shape_grad_component(shape_func,q_point,i)[j] +
		   fe_values.shape_grad_component(shape_func,q_point,j)[i]) / 2;
  return tmp;
}


// return vonmises stress
template <int dim>
inline
double
get_vonmises (const SymmetricTensor<2,dim> &stress)
{
  return sqrt( 0.5 * ( pow( (stress[0][0] - stress[1][1]), 2) 
		       + pow( (stress[1][1] - stress[2][2]), 2)
		       + pow( (stress[2][2] - stress[0][0]), 2)
		       + 6 * ( 
			      pow(stress[0][1], 2) 
			      + pow(stress[1][2], 2) 
			      + pow(stress[2][0], 2))
		       ) 
		 );
}


// return thermal strain i.e diagonal strain tensor
template <int dim>
inline
SymmetricTensor<2,dim>
  get_therm_strain (double &eps1, double &eps2, double &eps3)
{
  SymmetricTensor<2,dim> tmp;
  tmp[0][0] = eps1;
  tmp[1][1] = eps2;
  tmp[2][2] = eps3;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i+1; j<dim; ++j)
      tmp[i][j] = 0 ;
  return tmp;
}


// given the dispalcement gradients, return the strain tensor
template <int dim>
inline
SymmetricTensor<2,dim>
  get_strain (const std::vector<Tensor<1,dim> > &grad)
{
  Assert (grad.size() == dim, ExcInternalError());
  SymmetricTensor<2,dim> strain;
  for (unsigned int i=0; i<dim; ++i)
    strain[i][i] = grad[i][i];
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i+1; j<dim; ++j)
      strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
  return strain;
}


#endif
