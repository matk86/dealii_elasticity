#ifndef __ale_h
#define __ale_h

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <common.h>

namespace ale
{

using namespace dealii;


/*------------------------------------------------------------------
  class LameLambda
  ------------------------------------------------------------------*/
template <int dim>
class LameLambda :  public Function<dim>
{
public:
	LameLambda (double val);
	virtual double value (const Point<dim> &p,
			const unsigned int component = 0) const;
	virtual void value_list (const std::vector<Point<dim> > &points,
			std::vector<double> &value_list,
			const unsigned int component =0) const;
private:
	double ALE_LAMBDA;

};


/*------------------------------------------------------------------
  class LameMu
  ------------------------------------------------------------------*/
template <int dim>
class LameMu :  public Function<dim>
{
public:
	LameMu (double val);
	virtual double value (const Point<dim> &p,
			const unsigned int component = 0) const;
	virtual void value_list (const std::vector<Point<dim> > &points,
			std::vector<double>  &value_list,
			const unsigned int component = 0) const;
private:
	double ALE_MU;
};


/*------------------------------------------------------------------
  class BodyForce
  ------------------------------------------------------------------*/
template <int dim>
class BodyForce
{
public:
	BodyForce(Point<dim> p1, Point<dim> p2, Point<dim> p3,
		  std::vector<double> hkl, double line_tol, 
		  double fmag, double forceratio);
	bool isonline(const Point<dim> & p, const Point<dim> & pbegin,
			const Point<dim> & pend, double tol);
	std::vector<double> get_force(const Point<dim> & p);
	void force_value_list(const std::vector<Point<dim> > &points,
			std::vector<std::vector<double> >  &force_values);

private:
	double ALE_FORCE ;
	double ALE_FORCE_RATIO ;
	double ALE_TOL;
	Point<dim> ALE_P1;
	Point<dim> ALE_P2;
	Point<dim> ALE_P3;
	Tensor<1,3, double> ALE_T1;
	Tensor<1,3, double> ALE_T2;
	Tensor<1,3, double>  ALE_T3;
	Tensor<1,3, double> ALE_HKL;
};


/*------------------------------------------------------------------
  class ElasticProblem
  -------------------------------------------------------------------*/
template <int dim>
class ElasticProblem
{
public:
	ElasticProblem (std::string fname,
			double lambda_val, double mu_val, 
			Point<dim> p1, Point<dim> p2, Point<dim> p3,
			std::vector<double> hkl, 
			double line_tol, double fmag, double forceratio,
			std::vector<int> mesh_set, std::vector<int> refine_surf,
			unsigned int max_iter, double conv_tol, double precond_param,
			std::vector<int> dir_bc, std::vector<std::vector<bool> > dirbc_mask);
	~ElasticProblem ();
	void run ();
	void post_process();

private:
	void setup_system ();
	void setup_system_constr ();
	void assemble_and_solve ();
	void assemble_and_solve_constr ();
	void solve ();
	void refine_grid ();
	void output_results () const;

	Triangulation<dim>   triangulation;
	FESystem<dim>        fe;
	DoFHandler<dim>      dof_handler;
	MappingQ<dim>        mapping;
	double               elastic_energy;

	//ConstraintMatrix     hanging_node_constraints;
	ConstraintMatrix     constraints;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double>       solution;
	Vector<double>       strain_xx;
	Vector<double>       strain_zz;
	Vector<double>       system_rhs;

	std::string  ALE_MESH;
	int ALE_GLOBAL_REFINE;
	int ALE_LOCAL_REFINE ;
	int ALE_MAX_ITER ;
	double ALE_CONV_TOL ;
	double ALE_PRECOND_PARAM ;
	std::vector<int> ALE_DIR_BC;
	std::vector<int> ALE_REFINE_SURFACES;

	LameLambda<dim>  LLAMBDA;
	LameMu<dim>      LMU;
	BodyForce<dim>   LBDF;
	std::vector<std::vector<bool> > ALE_DIRBC_MASK;

	//static const SymmetricTensor<4,dim> stress_strain_tensor;
};


}

#endif
