std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals_R(int e2max, double hw, double Eclosure, double r12)
{
  IMSRGProfiler profiler;
  double t_start_pci = omp_get_wtime(); // profiling (s)
  std::unordered_map<uint64_t,double> IntList;
  int size=1000;
  std::cout<<"calculating integrals wrt dq..."<<std::endl;
  int maxn = e2max/2;
  int maxl = e2max;
  int maxnp = e2max/2;
  std::vector<uint64_t> KEYS;
  for (int S = 0; S<=1; S++)
  {
    for (int n=0; n<=maxn; n++)
    {
      for (int l=0; l<=maxl; l++)
      {
        int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,l) = IntHash(np,l,n,l), by construction
        for (int np=tempminnp; np<=maxnp; np++)
        {
          int minJ = abs(l-S);
          int tempmaxJ = l+S;
          for (int J = minJ; J<= tempmaxJ; J++)
          {
            uint64_t key = IntHash(n,l,np,l,S,J);
            KEYS.push_back(key);
            IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
          }
        }
      }
    }
  }

  gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(size);
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<KEYS.size(); i++)
  {
    uint64_t key = KEYS[i];
    uint64_t n,l,np,lp,S,J;
    IntUnHash(key, n,l,np,lp,S,J);
    IntList[key] = r12*r12*HO_Radial_psi(n, l, hw, r12)*HO_Radial_psi(np, lp, hw, r12)*integrate_dq_radial_GT(Eclosure,r12,size,t); // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
  }
  gsl_integration_glfixed_table_free(t);

  std::cout<<"...done calculating the integrals"<<std::endl;
  std::cout<<"IntList has "<<IntList.bucket_count()<<" buckets and a load factor "<<IntList.load_factor()
    <<", estimated storage ~= "<<((IntList.bucket_count() + IntList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<std::endl; // copied from (RS)
  profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci; // profiling (r)
  return IntList;
}