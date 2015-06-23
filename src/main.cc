#include <burger.h>
#include <ale.h>

int main (int argc, char* argv[])
{
  try
    {
      dealii::deallog.depth_console (0);
      
      std::string fname;
      double fmag, forceratio,line_tol, lambda_val, mu_val;
      dealii::Point<3> p1, p2, p3;
      std::vector<double> hkl(3);
      std::vector<int> mesh_refine, refine_surf;
      std::vector<int> dir_bc;
      std::vector<bool> which_comp1(3, true);
      std::vector<std::vector<bool> > which_comp;
      unsigned int max_iter;
      double conv_tol, precond_param;
      
      // input mesh
      if (argc < 2)
	{
	  std::cerr << "Usage: " << argv[0] << " <mesh file>" << std::endl;
	  return 1;
	}
      std::cout << "Mesh file : "<< argv[1] << std::endl;
      
      /*----------------------------------------------------------------
	ALE Parameters
	-----------------------------------------------------------------*/
      // elastic constant: lame parameters
      lambda_val = 30.0;//in GPa
      mu_val = 48.5; //in GPa
      
      // body force magnitude
      fmag = 10.0;
      forceratio = 10.0; // F_111/F_100
      
      //width of line along which the body force will be applied
      line_tol = 0.01;
      
      // points that define the plane
      p1[0] = 1.0;  p1[1]=1.0;  p1[2]=0.0;
      p2[0] = 0.0;  p2[1]=1.0;  p2[2]=1.0;
      p3[0] = 1.0;  p3[1]=0.0;  p3[2]=1.0;

      // normal to the plane
      hkl[0] = 1.0; hkl[1]=1.0; hkl[2] = 1.0;

      // global and local refinement steps
      mesh_refine.push_back(0); // global, 0 ==> no refinement
      mesh_refine.push_back(3); // local
      
      //surface indicators for local mesh refinement
      refine_surf.push_back(4);
      refine_surf.push_back(6);
      refine_surf.push_back(1);
      refine_surf.push_back(5);

      //boundary indicators of surfaces where dirichlet bc 
      //will be applied
      //note: zero displacement
      dir_bc.clear();
      which_comp.clear();
      // xy plane, z displacement set to 0
      dir_bc.push_back(2);
      which_comp1.clear();
      which_comp1.push_back(false);
      which_comp1.push_back(false);
      which_comp1.push_back(true);
      which_comp.push_back(which_comp1);
      // yz plane, x displacement set to 0
      dir_bc.push_back(7);
      which_comp1.clear();
      which_comp1.push_back(true);
      which_comp1.push_back(false);
      which_comp1.push_back(false);
      which_comp.push_back(which_comp1);
      // xz plane, y displacement set to 0
      dir_bc.push_back(3);
      which_comp1.clear();
      which_comp1.push_back(false);
      which_comp1.push_back(true);
      which_comp1.push_back(false);
      which_comp.push_back(which_comp1);

      // solver settings
      max_iter = 2000;
      conv_tol = 1e-6;
      precond_param = 0.6;
      /*----------------------------------------------------------------
	end Parameters
	-----------------------------------------------------------------*/
      
      ale::ElasticProblem<3> cuboct( argv[1],
				     lambda_val, mu_val,
				     p1, p2, p3,
				     hkl, line_tol, fmag, forceratio,
				     mesh_refine, refine_surf,
				     max_iter, conv_tol, precond_param,
				     dir_bc, which_comp);
      cuboct.run();

      // burger::ElasticProblem<3> nanoburger();
      // nanoburger.run()
      
      std::cout << "DONE ! "<< std::endl;
    }
  
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "--------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "--------------------------------------------------"
		<< std::endl;
      return 1;
    }
  
  catch (...)
    {
      std::cerr << std::endl << std::endl
		<< "--------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "--------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
