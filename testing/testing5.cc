  //Get an integral from the IntList cache or calculate it (parallelization dependent)
  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp,int J, double hw, std::string transition, double Eclosure, std::string src, std::unordered_map<uint64_t,double> &IntList)
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
    uint64_t key = IntHash(n,l,np,lp,J);
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
      std::unordered_map<uint64_t,double> AList = PrecalculateA(e2max,Eclosure,transition,size);
      integral = integrate_dq(n, l, np, lp,J,hw, transition, Eclosure, src,size, t,AList);
      gsl_integration_glfixed_table_free(t);
      if (omp_get_num_threads() >= 2)
      {
        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
        printf("   I shouldn't be here in GetIntegral(%d, %d, %d, %d, %d):   key =%llx   integral=%f\n",n,l,np,lp,J,key,integral);
        exit(EXIT_FAILURE);
      }
      IntList[key] = integral;
      
      return integral;
    }
  }