  //Get an integral from the IntList cache or calculate it (parallelization dependent)
  double GetM0nuIntegral_R(int e2max, int n, int l, int np, int lp, int S, int J, double hw, double Eclosure, double r12, std::unordered_map<uint64_t,double> &IntList)
  {
    int maxl = e2max;
    int maxnp = e2max/2;
    int maxlp = e2max;
    int order1 = decimalgen(n,l,np,lp,maxl,maxnp,maxlp);
    int order2 = decimalgen(np,lp,n,l,maxl,maxnp,maxlp); // notice I was careful here with the order of maxl,maxnp,maxlp to make proper comparison
    if (order1 > order2)
    {
      std::swap(n,np); // using symmetry IntHash(n,l,np,lp) = IntHash(np,lp,n,l)
      std::swap(l,lp); // " " " " "
    }
    // long int key = IntHash(n,l,np,lp); // if I ever get that version working...
    // std::cout<<"n ="<<n<<", l ="<<l<<", np = "<<np<<", S = "<<S<<", J = "<<J<<std::endl;
    uint64_t key = IntHash(n,l,np,lp,S,J);
    auto it = IntList.find(key);
    if (it != IntList.end()) // return what we've found
    {
      return it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      double integral;
      int size = 500;
      gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(size);
      integral = r12*r12*HO_Radial_psi(n, l, hw, r12)*HO_Radial_psi(np, lp, hw, r12)*integrate_dq_radial_GT(Eclosure,r12,size,t);
      gsl_integration_glfixed_table_free(t);
      if (omp_get_num_threads() >= 2)
      {
        std::cout << "DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!" << std::endl;
        std::cout << "   I shouldn't be here in GetIntegral(" << n << " , " << l << " , " << np << " , " << lp << " , " << S << " , " << J << "):   key = " << key << "   integral = " << integral << std::endl; 
//        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
//        printf("   I shouldn't be here in GetIntegral(%d, %d, %d, %d, %d, %d):   key =%lx   integral=%f\n",n,l,np,lp,S,J,key,integral);
        exit(EXIT_FAILURE);
      }
      IntList[key] = integral;

      return integral;
    }
  }