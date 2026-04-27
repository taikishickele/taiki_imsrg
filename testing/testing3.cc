std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals(int e2max, double hw, std::string transition, double Eclosure, std::string src)
  {
    IMSRGProfiler profiler;
    double t_start_pci = omp_get_wtime(); // profiling (s)
    std::unordered_map<uint64_t,double> IntList;
    std::unordered_map<uint64_t,double> AList;
    int size=500;
    std::cout<<"calculating integrals wrt dp and dpp..."<<std::endl;
    if (transition != "C")
    {
      AList = PrecalculateA(e2max,Eclosure,transition,size);
      std::cout<<"Done precomputing A's."<<std::endl;
    }
    int maxn = e2max/2;
    int maxl = e2max;
    int maxnp = e2max/2;
    std::vector<uint64_t> KEYS;
    if (transition == "F" or transition == "GT")
    {
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
                uint64_t key = IntHash(n,l,np,l,J);
                KEYS.push_back(key);
                IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
              }
            }
          }
        }
      }
      
    }
    else if (transition == "T")
    {
      for (int n=0; n<=maxn; n++)
      {
        for (int l=1; l<=maxl; l++)
        {
          int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
          //int tempminnp = 0;
          for (int np=tempminnp; np<=maxnp; np++)
          {
            int tempminlp = (n==np ? l : 1); // NOTE: need not start from 'int lp=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
            int maxlp = std::min(l+2,maxl);
            for (int lp = tempminlp; lp<=maxlp; lp++)
            { 
              if ((abs(lp-l) != 2) and (abs(lp-l) != 0)) continue;
              int minJ = std::max(abs(l-1),abs(lp-1));
              int tempmaxJ = std::min(l+1,lp+1);
              for (int J = minJ; J<= tempmaxJ; J++)
              {
                uint64_t key = IntHash(n,l,np,lp,J);
                KEYS.push_back(key);
                IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
              }
            }
          }
        }
      }
    }
    else if (transition =="C")
    {
      int S = 0;
      for (int n=0; n<=maxn; n++)
      {
        int l = 0;
        int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,l) = IntHash(np,l,n,l), by construction
        for (int np=tempminnp; np<=maxnp; np++)
        {
          int J = 0;
          uint64_t key = IntHash(n,l,np,l,J);
          KEYS.push_back(key);
          IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
        }
      }
    }

    gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(size);
    #pragma omp parallel for schedule(dynamic, 1)// this works as long as the gsl_function handle is within this for-loop
    for (size_t i=0; i<KEYS.size(); i++)
    {
      uint64_t key = KEYS[i];
      uint64_t n,l,np,lp,J;
      IntUnHash(key, n,l,np,lp,J);
      IntList[key] = integrate_dq(n,l,np,lp,J,hw,transition,Eclosure,src,size,t,AList); // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
      // Uncomment if you want to verify integrals values
      // std::stringstream intvalue;
      // intvalue<<n<<" "<<l<<" "<<np<<" "<<lp<<" "<<J<<" "<<IntList[key]<<std::endl;
      // std::cout<<intvalue.str();
      
    }
    gsl_integration_glfixed_table_free(t);
    
    std::cout<<"...done calculating the integrals"<<std::endl;
    std::cout<<"IntList has "<<IntList.bucket_count()<<" buckets and a load factor "<<IntList.load_factor()
      <<", estimated storage ~= "<<((IntList.bucket_count() + IntList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<std::endl; // copied from (RS)
    profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci; // profiling (r)
    return IntList;
  }