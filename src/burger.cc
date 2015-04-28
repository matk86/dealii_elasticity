#include <burger.h>
#include <utility.h>

namespace burger
{
  using namespace dealii;

  // ElasticProblem::ElasticProblem
  template <int dim>
  ElasticProblem<dim>::ElasticProblem () : fe (FE_Q<dim>(1), dim), 
					   dof_handler(triangulation),
					   mapping(3), 
					   elastic_energy(0)
  {}

  //ElasticProblem::~ElasticProblem
  template <int dim>
  ElasticProblem<dim>::~ElasticProblem ()
  {
    dof_handler.clear ();
  }

  //ElasticProblem::setup_system
  template <int dim>
  void ElasticProblem<dim>::setup_system ()
  {
    //allocate spce for the dofs with the chosen lagrange 
    //shpae functions
    dof_handler.distribute_dofs (fe);

    //sparsity and hanging nodes
    //hanging_node_constraints.clear ();
    //    DoFTools::make_hanging_node_constraints (dof_handler, 
    //					     hanging_node_constraints);
    //hanging_node_constraints.close ();
    //sparsity_pattern.reinit (dof_handler.n_dofs(), 
    //			     dof_handler.n_dofs(), 
    //			     dof_handler.max_couplings_between_dofs());
    //DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    //hanging_node_constraints.condense (sparsity_pattern);
    //sparsity_pattern.compress();

    //initialize
    //system_matrix.reinit (sparsity_pattern);
    //CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    //DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    //sparsity_pattern.copy_from(c_sparsity);
    DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    sparsity_pattern.copy_from(c_sparsity);


    system_matrix.reinit (sparsity_pattern);
    solution.reinit (dof_handler.n_dofs());
    strain_xx.reinit (triangulation.n_active_cells());
    strain_zz.reinit (triangulation.n_active_cells());

    system_rhs.reinit (dof_handler.n_dofs());
  }

  //ElasticProblem::assemble_and_solve
  template <int dim>
  void ElasticProblem<dim>::assemble_and_solve ()
  {
    //setup the gaussian quatrature formula
    //quadrature rule should have at least the order of the 
    //boundary approximation(i.e the mapping)
    const unsigned int gauss_degree = 
      std::max (
		static_cast<unsigned int>(std::ceil(1.*(mapping.get_degree()+1)/2)), 2U);
    QGauss<dim>  quadrature_formula(gauss_degree);//2);

    //setup the shape func evaluator 
    FEValues<dim> fe_values (mapping, 
			     fe, 
			     quadrature_formula,
			     update_values   | update_gradients |
			     update_quadrature_points | 
			     update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    LameLambda<dim>         lambda;
    LameMu<dim>             mu;
    Temperature<dim>        tempe;
    // ConstantFunction<dim> lambda(1.), mu(1.);
    std::vector<double>     lambda_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);
    std::vector<double>     t_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = 
      dof_handler.begin_active(), endc = dof_handler.end();

    //setup the cell matrix and the cell rhs for each cell
    //loop over the cells and the cell dofs
    for (; cell!=endc; ++cell)
      {
	cell_matrix = 0;
        cell_rhs = 0;

	fe_values.reinit (cell);

	lambda.value_list(fe_values.get_quadrature_points(), 
			  lambda_values);
	mu.value_list(fe_values.get_quadrature_points(), 
		      mu_values);
	tempe.value_list(fe_values.get_quadrature_points(), 
			 t_values);

        //cell matrix
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
		      +=
		      (
			(fe_values.shape_grad(i,q_point)[component_i] *
			 fe_values.shape_grad(j,q_point)[component_j] *
			 lambda_values[q_point])
			+
			(fe_values.shape_grad(i,q_point)[component_j] *
			 fe_values.shape_grad(j,q_point)[component_i] *
			 mu_values[q_point])
			+
			((component_i == component_j) ?
			 (fe_values.shape_grad(i,q_point) *
			  fe_values.shape_grad(j,q_point) *
			  mu_values[q_point])  :
			 0)
		      )
		      *
		      fe_values.JxW(q_point);
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
	        //cell_rhs(i) += fe_values.shape_value(i,q_point) *
	        //	     rhs_values[q_point](component_i) *
	        //	     fe_values.JxW(q_point) + 
	        //           10.0 * 
		//      fe_values.shape_grad(i,q_point)[component_i] * 
		//      fe_values.JxW(q_point);

                //zero body force, only thermal contribution, 
		//the coeffecient = (3*lamda + 2* mu)*alpha*delta T
                /*if (component_i == 2)
		    cell_rhs(i) += 2 * lambda_values[q_point] * 
                                   t_values[q_point] * 
                                   fe_values.shape_grad(i,q_point)[component_i] * 
                                   fe_values.JxW(q_point);
	        else
	          cell_rhs(i) += 2 * (lambda_values[q_point] + mu_values[q_point]) * 
                                 t_values[q_point] * 
                                 fe_values.shape_grad(i,q_point)[component_i] * 
                                 fe_values.JxW(q_point);*/
		 cell_rhs(i) += 
		   (2*lambda_values[q_point] + 
		    2*mu_values[q_point]) * 
		   t_values[q_point] * 
		   fe_values.shape_grad(i,q_point)[component_i] * 
		   fe_values.JxW(q_point);
	      }
	  }

        //assemble the global matrix and rhs
	cell->get_dof_indices (local_dof_indices);

	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      system_matrix.add(local_dof_indices[i], 
				local_dof_indices[j], 
				cell_matrix(i,j));
              system_rhs(local_dof_indices[i]) += cell_rhs(i);
	  }
      }//end of the loop over cell

    //hanging_node_constraints.condense (system_matrix);
    //hanging_node_constraints.condense (system_rhs);

    //Apply essential boundary conditions
    //zero displacemnt on the surface with boundary indicator = 1
    std::map<unsigned int,double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
					     1,
					     ZeroFunction<dim>(dim),
					     boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, 
				       system_matrix, 
				       solution, system_rhs);

    solve();

    //PostProcess
    std::vector<std::vector<Tensor<1,dim> > > 
      solution_grads(quadrature_formula.size(),
		     std::vector<Tensor<1,dim> >(dim));

    double ztop, xlt=-0.83, xrt=0.83;
    Vector<double> disp(3);
    for(ztop = 0.54; ztop<=0.55; ztop += 0.01)
      {
        VectorTools::point_value(dof_handler, 
				 solution,Point<3>(0., 0., ztop), 
				 disp);
	std::cout << "soln at " << 0.0 <<" " << 0.0 << " " 
		  <<ztop << " " <<  disp(0) << " "<< disp(1) 
		  << " "<<disp(2) << " "<< "\n";
        VectorTools::point_value(dof_handler, 
				 solution,Point<3>(xlt, 0., ztop), 
				 disp);
	std::cout << "soln at " << xlt <<" " << 0.0 << " " 
		  <<ztop << " " <<  disp(0) << " "<< disp(1) 
		  << " "<<disp(2) << " "<< "\n";
        VectorTools::point_value(dof_handler, 
				 solution,Point<3>(xrt, 0., ztop), 
				 disp);
	std::cout << "soln at " << xrt <<" " << 0.0 << " " <<ztop 
		  << " " <<  disp(0) << " "<< disp(1) << " "<<disp(2)
		  << " "<< "\n";
      }

    elastic_energy = 0.0;

    unsigned int i = 0;

    for (typename DoFHandler<dim>::active_cell_iterator cell = 
	   dof_handler.begin_active(); 
	 cell != dof_handler.end(); 
	 ++cell)
      {
        //std::cout <<"face_dofhandler" << cell->face(0)->center()[0];
        fe_values.reinit (cell);

	lambda.value_list(fe_values.get_quadrature_points(), 
			  lambda_values);
	mu.value_list(fe_values.get_quadrature_points(), mu_values);
	tempe.value_list(fe_values.get_quadrature_points(), t_values);
        fe_values.get_function_gradients (solution,solution_grads);

        double elastic_energy_cell = 0.0;
        double xx_sum = 0.0;
        double zz_sum = 0.0;

        for (unsigned int q=0; q<quadrature_formula.size(); ++q)
          {
            const SymmetricTensor<2,dim> eps_ij = 
	      get_strain (solution_grads[q]);

            double eps1 = t_values[q];
            SymmetricTensor<2,dim> eps_therm_ij;
	    eps_therm_ij[0][0] = 
	      eps1; eps_therm_ij[1][1] = 
		      eps1; eps_therm_ij[2][2] = eps1;
            for (unsigned int ii=0; ii<dim; ++ii)
              for (unsigned int jj=ii+1; jj<dim; ++jj)
                eps_therm_ij[ii][jj] = 0 ;
              //const SymmetricTensor<2,dim> eps_therm_ij = 
	      //get_therm_strain (eps1,eps1,eps1);

              const SymmetricTensor<2,dim> eps_elast_ij = 
		eps_ij - eps_therm_ij;

	      xx_sum += eps_elast_ij[0][0];
	      zz_sum += eps_elast_ij[2][2];

              const SymmetricTensor<4,dim> stress_strain_tensor = 
		get_stress_strain_tensor<dim> (lambda_values[q], 
					       mu_values[q]);

              const SymmetricTensor<2,dim> stress =
		stress_strain_tensor * eps_elast_ij;

              elastic_energy_cell +=  
		eps_elast_ij * 
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
    SolverControl           solver_control (2000, 1e-12);
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

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
  ElasticProblem<dim>::output_results(const unsigned int cycle) const
  {
    std::string filename = "solution-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

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
    for (unsigned int cycle=0; cycle<1; ++cycle)
      {
	std::cout << "Cycle " << cycle << ':' << std::endl;

        //intialize the geometry and the mesh
	if (cycle == 0)
	  {
	    GridIn<dim> grid_in;
	    GridOut grid_out;
	    grid_in.attach_triangulation (triangulation);
	    std::ifstream input_file("cuboct.msh");
	    std::ofstream output_file("cuboct_mesh.vtk");
	    std::ofstream output_file_msh("cuboct_mesh.msh");
	    grid_in.read_msh(input_file);
	    triangulation.refine_global(3);
	    grid_out.write_vtk(triangulation,output_file);
	    grid_out.write_msh(triangulation,output_file_msh);
	    std::cout <<"done";
	    exit(0);

	    //GridGenerator::hyper_ball (triangulation, 
	    //Point<dim>(0.0,0.0,0.0), 1.0);
            static const HyperBallBoundary<dim> boundary;

            triangulation.set_boundary (0, boundary);

      	    triangulation.refine_global(5);

	    for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin(); cell != triangulation.end(); ++cell)
	      {
	        //std::cout <<"face" << cell->face(0)->center()[0] 
		//<< "\n";//<< "  "<<cell->face(f)->center()[1] 
		//<< "  "<<cell->face(f)->center()[2] <<"\n" ;
	        for (unsigned int f=0; 
		     f<GeometryInfo<dim>::faces_per_cell; 
		     ++f)
		  {
		    if (cell->face(f)->at_boundary())
		      {
		        //if ((cell->face(f)->center()[0] > 0.9) || 
			//(cell->face(f)->center()[1] > 0.9))
		        //if ((cell->face(f)->center()[2] > -0.1) && 
			//(cell->face(f)->center()[2] < 0.1))
		        //if (fabs(cell->face(f)->center()[2]) <=1e-6 )
                       if ( (cell->face(f)->center()[dim-1] + 1.)<1e-6)
			 // || (cell->face(f)->center()[dim-1] > 0.95))
			  {
               		    cell->face(f)->set_all_boundary_indicators(1);
			  }
		      }
		  }	
	      }

            //refine the region around the interface
	    for (unsigned int step=0; step<2; ++step)
	      {
	        typename Triangulation<dim>::active_cell_iterator 
		  acell= triangulation.begin_active(), 
		  endc = triangulation.end();

	        for (; acell!=endc; ++acell)
	          for (unsigned int v=0; 
		       v < GeometryInfo<dim>::vertices_per_cell; ++v)
	   	    {
                      //if (std::fabs(acell->vertex(v)[dim-1]) < 0.2)
                      if ( ((acell->vertex(v)[dim-1] > 0.45) && 
			    (acell->vertex(v)[dim-1] < 0.65)) || 
			   ((acell->vertex(v)[dim-1] < -0.45) && 
			    (acell->vertex(v)[dim-1] > -0.65)))
		      {
		        acell->set_refine_flag ();
		        break;
		      }
		    } 
	        triangulation.execute_coarsening_and_refinement ();
	      }
	 }
	else
      	  //triangulation.refine_global (1);
	  refine_grid ();

	std::cout << "   Number of active cells: "
		  << triangulation.n_active_cells() << std::endl;

	setup_system ();

	std::cout << "   Number of degrees of freedom: " 
		  << dof_handler.n_dofs() << std::endl;

	assemble_and_solve ();
	//solve ();
        //post_process();
	output_results (cycle);
      }
  }



/*----------------------------------------------------------------------------
  ElasticProblem::post_process
----------------------------------------------------------------------------*/
/*  template <int dim>
  void ElasticProblem<dim>::post_process ()
  {
    const unsigned int gauss_degree = std::max (static_cast<unsigned int>(std::ceil(1.*(mapping.get_degree()+1)/2)),2U);
    QGauss<dim>  quadrature_formula(gauss_degree);//2);

    const unsigned int   n_q_points    = quadrature_formula.size();

    // ConstantFunction<dim> lambda(1.), mu(1.);
    LameLambda<dim> lambda;
    LameMu<dim>  mu;
    std::vector<double>     lambda_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    FEValues<dim> fe_values (fe, quadrature_formula,update_values | update_gradients | update_JxW_values);

    std::vector<std::vector<Tensor<1,dim> > > solution_grads (quadrature_formula.size(), std::vector<Tensor<1,dim> >(dim));
    //std::vector< Vector<double > > sol_bou;

    elastic_energy = 0.0;
    unsigned int i = 0;

    for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
      {
        //std::cout <<"face_dofhandler" << cell->face(0)->center()[0];
        fe_values.reinit (cell);

        lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
        mu.value_list     (fe_values.get_quadrature_points(), mu_values);

        fe_values.get_function_grads (solution,solution_grads);

        double elastic_energy_cell = 0.0;
        double trace_sum = 0.0;

        for (unsigned int q=0; q<quadrature_formula.size(); ++q)
          {
            const SymmetricTensor<2,dim> eps_ij = get_strain (solution_grads[q]);
	    trace_sum += trace(eps_ij);

            const SymmetricTensor<4,dim> stress_strain_tensor = get_stress_strain_tensor<dim> (lambda_values[q], mu_values[q]);

            const SymmetricTensor<2,dim> stress = stress_strain_tensor * eps_ij;
            elastic_energy_cell +=  eps_ij * stress * fe_values.JxW(q);
	  }
	  strain_trace(i) = trace_sum/quadrature_formula.size();
	  ++i;    
	  elastic_energy += elastic_energy_cell;
       }
    std::cout << "  elastic energy = " << elastic_energy << "\n";
    }*/

  // explicit initialization of the template
  template class ElasticProblem<3>;

}
