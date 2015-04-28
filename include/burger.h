#ifndef __burger_h
#define __burger_h

#include <deal.II/base/function.h>
#include <common.h>

namespace burger
{
  using namespace dealii;

  /*------------------------------------------------------------------
  class LameLambda
  ------------------------------------------------------------------*/
  template <int dim>
  class LameLambda :  public Function<dim>
  {
    public:
      LameLambda ();
      virtual double value (const Point<dim> &p, 
			    const unsigned int component = 0) const;
      virtual void value_list (const std::vector<Point<dim> > &points,
			       std::vector<double> &value_list, 
			       const unsigned int component =0) const;
  };


  /*------------------------------------------------------------------
  class LameMu
  ------------------------------------------------------------------*/
  template <int dim>
  class LameMu :  public Function<dim>
  {
    public:
    LameMu ();
    virtual double value (const Point<dim> &p, 
			  const unsigned int component = 0) const;
    virtual void value_list (const std::vector<Point<dim> > &points, 
			     std::vector<double>  &value_list, 
			     const unsigned int component = 0) const;
  };


  /*------------------------------------------------------------------
  class Temperature
  ------------------------------------------------------------------*/

  template <int dim>
  class Temperature :  public Function<dim>
  {
  public:
    Temperature ();
    virtual double value (const Point<dim> &p, 
			  const unsigned int component = 0) const;
    virtual void value_list (const std::vector<Point<dim> > &points, 
			     std::vector<double>   &value_list, 
			     const unsigned int component = 0) const;
  };


  /*------------------------------------------------------------------
  class ElasticProblem
  -------------------------------------------------------------------*/
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem ();
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
    void output_results (const unsigned int cycle) const;
    
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

    //static const SymmetricTensor<4,dim> stress_strain_tensor;
  };


}

#endif
