#include <ale.h>
#include <utility.h>

namespace ale
{
using namespace dealii;

// ElasticProblem::ElasticProblem
// intialize the FESystem  with first order lagrange polynomial
// shape function
// third order mapping : needed for mapping curvilinear boundaries
template <int dim>
ElasticProblem<dim>::ElasticProblem (std::string fname, 
				     double lambda_val, double mu_val,
				     Point<dim> p1, Point<dim> p2, Point<dim> p3,
				     std::vector<double> hkl,
				     double line_tol, double fmag, double forceratio,
				     std::vector<int> mesh_set, 
				     std::vector<int> refine_surf,
				     unsigned int max_iter, double conv_tol, 
				     double precond_param,
				     std::vector<int> dir_bc, 
				     std::vector<std::vector<bool> > dirbc_mask)
  :
  fe (FE_Q<dim>(1), dim),
  dof_handler(triangulation),
  mapping(3),
  elastic_energy(0),
  LLAMBDA(lambda_val),
  LMU(mu_val),
  LBDF(p1, p2, p3, hkl, line_tol, fmag, forceratio)
{
  ALE_MESH = fname;
  ALE_GLOBAL_REFINE = mesh_set[0];
  ALE_LOCAL_REFINE = mesh_set[1];
  ALE_REFINE_SURFACES =  refine_surf;
  ALE_MAX_ITER = max_iter;
  ALE_CONV_TOL = conv_tol;
  ALE_PRECOND_PARAM = precond_param;
  ALE_DIR_BC =  dir_bc;
  ALE_DIRBC_MASK = dirbc_mask;
}

// ElasticProblem::~ElasticProblem
template <int dim>
ElasticProblem<dim>::~ElasticProblem ()
{
	dof_handler.clear ();
}

// ElasticProblem::setup_system
template <int dim>
void
ElasticProblem<dim>::setup_system ()
{
	//allocate space for the dofs with the chosen FEsystem object
	// FEsystem object describes the number of unknowns on each
	// node as well the degree of shape functions
	dof_handler.distribute_dofs (fe);

	//sparsity and hanging nodes
	//hanging_node_constraints.clear ();
	//DoFTools::make_hanging_node_constraints (dof_handler,
	//				     hanging_node_constraints);
	//hanging_node_constraints.close ();
	//sparsity_pattern.reinit (dof_handler.n_dofs(),
	//			     dof_handler.n_dofs(),
	//			     dof_handler.max_couplings_between_dofs());
	//DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
	//hanging_node_constraints.condense (sparsity_pattern);
	//sparsity_pattern.compress();
	DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
	sparsity_pattern.copy_from(c_sparsity);

	system_matrix.reinit (sparsity_pattern);
	system_rhs.reinit (dof_handler.n_dofs());
	solution.reinit (dof_handler.n_dofs());
	strain_xx.reinit (triangulation.n_active_cells());
	strain_zz.reinit (triangulation.n_active_cells());
}

//ElasticProblem::assemble_and_solve
template <int dim>
void
ElasticProblem<dim>::assemble_and_solve ()
{
	//setup the gaussian quadrature formula
	//quadrature rule should have at least the order of the
	//boundary approximation(i.e the mapping)
	const unsigned int gauss_degree =
			std::max (
				  static_cast<unsigned int>(std::ceil(1.*(mapping.get_degree()+1)/2)), 2U);
	QGauss<dim>  quadrature_formula(gauss_degree);//2);

	// setup the shape function evaluator
	// evaluates all shape function related values at the
	// quadrature points
	FEValues<dim> fe_values (
			mapping,
			fe,
			quadrature_formula,
			update_values   | update_gradients |
			update_quadrature_points |
			update_JxW_values
	);

	const unsigned int  dofs_per_cell = fe.dofs_per_cell;
	// number of quadrature points in the cell
	const unsigned int  n_q_points    = quadrature_formula.size();

	FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs (dofs_per_cell);

	std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	// ConstantFunction<dim> lambda(1.), mu(1.);
	std::vector<double>     lambda_values (n_q_points);
	std::vector<double>     mu_values (n_q_points);
	std::vector<std::vector<double> >  force_values (n_q_points, std::vector<double>(3,0.0));

	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();

	//setup the cell matrix and the cell rhs for each cell
	//loop over the cells and the cell dofs
	for (; cell!=endc; ++cell)
	{
		cell_matrix = 0;
		cell_rhs = 0;

		fe_values.reinit (cell);

		this->LLAMBDA.value_list(fe_values.get_quadrature_points(),
				lambda_values);
		this->LMU.value_list(fe_values.get_quadrature_points(),
				mu_values);
		force_values.clear();
		this->LBDF.force_value_list(fe_values.get_quadrature_points(),
				force_values);
		//exit(0);

		//cell matrix
		// refer step:8 for linear elasticity stiffness matrix and
		// the rhs vector
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			const unsigned int component_i =
					fe.system_to_component_index(i).first;

			for (unsigned int j=0; j<dofs_per_cell; ++j)
			{
				const unsigned int component_j =
						fe.system_to_component_index(j).first;

				for (unsigned int q_point=0; q_point<n_q_points;
						++q_point)
				{
					cell_matrix(i,j)
					  += ( ( fe_values.shape_grad(i,q_point)[component_i] *
						 fe_values.shape_grad(j,q_point)[component_j] *
						 lambda_values[q_point] )
					       +
					       ( fe_values.shape_grad(i,q_point)[component_j] *
						 fe_values.shape_grad(j,q_point)[component_i] *
						 mu_values[q_point])
					       +
					       ( (component_i == component_j) ?
						 (fe_values.shape_grad(i,q_point) *
						  fe_values.shape_grad(j,q_point) *
						  mu_values[q_point])
						 :
						 0 )
					       ) * fe_values.JxW(q_point);
				}
			}
		}

		//rhs 
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			const unsigned int component_i =
					fe.system_to_component_index(i).first;

			for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
			{
			  cell_rhs(i) +=
			    force_values[q_point][component_i] *
			    fe_values.shape_value_component(i,q_point, component_i) *
			    fe_values.JxW(q_point);
			}
		}

		//assemble the global matrix and rhs
		cell->get_dof_indices (local_dof_indices);

		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  system_matrix.add( local_dof_indices[i],
					     local_dof_indices[j],
					     cell_matrix(i,j) );
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}//end of the loop over cell

	//hanging_node_constraints.condense (system_matrix);
	//hanging_node_constraints.condense (system_rhs);

	//Apply essential boundary conditions
	//zero displacement on the surfaces tagged by the  boundary indicators
	std::map<unsigned int,double> boundary_values;
	int it_count = 0;
	for(std::vector<int>::iterator it = this->ALE_DIR_BC.begin(); 
			it != this->ALE_DIR_BC.end();
			++it) {
	  std::cout << "Applying Dirichlet bc(zero displacement) on surface : " 
		    << *it
		    << " in directions, 0 : " << this->ALE_DIRBC_MASK[it_count][0] 
		    << " 1 : " << this->ALE_DIRBC_MASK[it_count][1]
		    << " 2 : " << this->ALE_DIRBC_MASK[it_count][2]
		    << std::endl ;
	  VectorTools::interpolate_boundary_values(
						   dof_handler,
						   *it,
						   ZeroFunction<dim>(dim),
						   boundary_values,
						   ComponentMask(this->ALE_DIRBC_MASK[it_count]));
	  MatrixTools::apply_boundary_values(
					     boundary_values,
					     system_matrix,
					     solution,
					     system_rhs );
	  ++it_count;
	}

	//solve
	std::cout << "Solving ...  " << std::endl;
	solve();

	//PostProcess
	// compute xx, zz strain as well as the elastic energy
	std::cout << "Post-processing ...  " << std::endl;
	std::vector<std::vector<Tensor<1,dim> > > solution_grads(
			quadrature_formula.size(),
			std::vector<Tensor<1,dim> >(dim) );

	elastic_energy = 0.0;

	unsigned int i = 0;

	for (typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active();
			cell != dof_handler.end();
			++cell)
	{
		//std::cout <<"face_dofhandler" << cell->face(0)->center()[0];
		fe_values.reinit (cell);

		this->LLAMBDA.value_list(fe_values.get_quadrature_points(),
				lambda_values);
		this->LMU.value_list(fe_values.get_quadrature_points(), mu_values);
		fe_values.get_function_gradients (solution,solution_grads);

		double elastic_energy_cell = 0.0;
		double xx_sum = 0.0;
		double zz_sum = 0.0;

		for (unsigned int q=0; q<quadrature_formula.size(); ++q)
		{
			const SymmetricTensor<2,dim> eps_ij =
					get_strain (solution_grads[q]);

			xx_sum += eps_ij[0][0];
			zz_sum += eps_ij[2][2];

			const SymmetricTensor<4,dim> stress_strain_tensor =
					get_stress_strain_tensor<dim> (lambda_values[q],
							mu_values[q]);

			const SymmetricTensor<2,dim> stress =
					stress_strain_tensor * eps_ij;

			elastic_energy_cell +=
					eps_ij *
					stress *
					fe_values.JxW(q);
		}

		strain_xx(i) = xx_sum / quadrature_formula.size();
		strain_zz(i) = zz_sum / quadrature_formula.size();

		++i;
		elastic_energy += elastic_energy_cell;
	}
	std::cout << "  elastic energy = " << elastic_energy << "\n";
}


// ElasticProblem::solve
// preconditioned CG
template <int dim>
void
ElasticProblem<dim>::solve ()
{
	SolverControl  solver_control (this->ALE_MAX_ITER, this->ALE_CONV_TOL);
	SolverCG<>   cg (solver_control);

	//PreconditionSSOR<> preconditioner;
	//preconditioner.initialize(system_matrix, 1.2);
	PreconditionJacobi<> preconditioner;
	preconditioner.initialize(system_matrix, this->ALE_PRECOND_PARAM);


	cg.solve (system_matrix, solution, system_rhs, preconditioner);
	//hanging_node_constraints.distribute (solution);
}

// ElasticProblem::refine_grid
// local mesh refinement based on error estimation
template <int dim>
void
ElasticProblem<dim>::refine_grid ()
{
	Vector<float>
	estimated_error_per_cell(triangulation.n_active_cells());

	typename FunctionMap<dim>::type neumann_boundary;
	KellyErrorEstimator<dim>::estimate (dof_handler,
			QGauss<dim-1>(2),
			neumann_boundary,
			solution,
			estimated_error_per_cell);

	GridRefinement::refine_and_coarsen_fixed_number(triangulation,
			estimated_error_per_cell,
			0.3, 0.03);

	triangulation.execute_coarsening_and_refinement ();
}

// ElasticProblem::output_results
// vtk output
template <int dim>
void
ElasticProblem<dim>::output_results() const
{
	std::string filename = "solution";
	filename += ".vtk";
	std::ofstream output (filename.c_str());

	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);

	std::vector<std::string> solution_names1(dim,"displacement");
	std::vector<std::string> solution_names2, solution_names3;;
	solution_names2.push_back ("strain_xx");
	solution_names3.push_back ("strain_zz");

	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	data_component_interpretation(dim,
			DataComponentInterpretation::component_is_part_of_vector);

	data_out.add_data_vector (solution,
			solution_names1,
			DataOut<dim>::type_dof_data,
			data_component_interpretation);
	data_out.add_data_vector (strain_xx,
			solution_names2,
			DataOut<dim>::type_cell_data);
	data_out.add_data_vector (strain_zz,
			solution_names3,
			DataOut<dim>::type_cell_data);
	data_out.build_patches ();

	data_out.write_vtk (output);
}

// ElasticProblem::run
template <int dim>
void
ElasticProblem<dim>::run ()
{
	GridIn<dim> grid_in;
	GridOut grid_out;
	std::string mesh_file_vtk;

	mesh_file_vtk = "mesh.vtk";

	grid_in.attach_triangulation (triangulation);

	std::cout << "reading in ... " << this->ALE_MESH << std::endl;
	std::ifstream input_file(this->ALE_MESH.c_str());
	if (!input_file) {
		std::cout << "Error while opening the file" << std::endl;
	}
	grid_in.read_msh(input_file);

	std::cout << "Number of active cells: "
			<< triangulation.n_active_cells() << std::endl;

	//refine all elements once
	if (this->ALE_GLOBAL_REFINE != 0)
	  {
	    std::cout << "Global refinement ...  " << this->ALE_GLOBAL_REFINE 
		      << " times" << std::endl;
	    triangulation.refine_global(this->ALE_GLOBAL_REFINE);
	    std::cout << "Number of active cells: "
		      << triangulation.n_active_cells() << std::endl;
	  }
	else
	    std::cout << "NO global refinement " << std::endl;
	//std::cout << "writing the mesh in vtk format ... "
	//			<< mesh_file_vtk << std::endl;
	//std::ofstream output_file(mesh_file_vtk.c_str());
	//grid_out.write_vtk(triangulation, output_file);
	//std::ofstream output_file_msh("mesh.msh");
	//grid_out.write_msh(triangulation,output_file_msh);

	std::cout << "Local refinement ...  " << std::endl;
	for (std::vector<int>::iterator s_it = this->ALE_REFINE_SURFACES.begin(); 
	     s_it != this->ALE_REFINE_SURFACES.end(); ++s_it)
	  std::cout << "Plane : " << *s_it << std::endl;

	std::vector<int>::iterator surface_it;
	//refine 3 times 111 plane and the plane around it
	for (unsigned int step=0; step<this->ALE_LOCAL_REFINE; ++step)
	{
		typename Triangulation<dim>::active_cell_iterator
		cell= triangulation.begin_active(), 
		endc = triangulation.end();

		for (; cell!=endc; ++cell)
			for (unsigned int f=0;
					f<GeometryInfo<dim>::faces_per_cell;
					++f)
			{
				if (cell->face(f)->at_boundary())
				{
					surface_it = std::find (this->ALE_REFINE_SURFACES.begin(),
							this->ALE_REFINE_SURFACES.end(),
							(int)cell->face(f)->boundary_indicator());
					if (surface_it != this->ALE_REFINE_SURFACES.end())
					{
						cell->set_refine_flag ();
						break;
						//exit(0);
					}
				}
			}
		//exit(0);
		triangulation.execute_coarsening_and_refinement ();
	}

	//write the refined mesh in vtk format
	//mesh_file_vtk = "mesh_refined.vtk";
	//std::cout << "writing the refined mesh in vtk format ... "
	//		<< mesh_file_vtk << std::endl;
	//std::ofstream output_file_refined(mesh_file_vtk.c_str());
	//grid_out.write_vtk(triangulation, output_file_refined);

	std::cout << "Number of active cells: "
			<< triangulation.n_active_cells() << std::endl;
	//exit(0);
	// initialize the fe dofs as well as the sytem matrix,
	// solution vector and the system rhs
	std::cout << "Setting up the system ...  " << std::endl;
	setup_system ();

	std::cout << "Number of degrees of freedom: "
			<< dof_handler.n_dofs() << std::endl;
	//exit(0);

	std::cout << "Assembling ...  " << std::endl;
	assemble_and_solve ();

	//post_process();
	std::cout << "Output to vtk files ... " << std::endl;
	output_results ();
}

// explicit initialization of the template
template class ElasticProblem<3>;

}
/*
					if ( (cell->face(f)->boundary_indicator() == 4) ||
							(cell->face(f)->boundary_indicator() == 6) ||
							(cell->face(f)->boundary_indicator() == 1) ||
							(cell->face(f)->boundary_indicator() == 5) )
 */
