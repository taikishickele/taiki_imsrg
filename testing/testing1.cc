Operator GamowTeller_old(ModelSpace& modelspace, double Eclosure, std::string src) // The old implementation of the Gamow-Teller operator w/o the PWD class
  {
    bool reduced = true;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    std::string transition = "GT";
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    M0nuGT_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [fm]
    const double prefact = Rnuc/(PI*PI); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [fm]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency
    M0nuGT_TBME.profiler.timer["M0nuGT_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout<<"calculating M0nu TBMEs..."<<std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto& itmat : M0nuGT_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      int J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat; // set below based on "reduced" variable
      if (reduced == false)
      {
        Jhat = 1.0; // for non-reduced elements, to compare with JE
      }
      else //if (reduced == "true")
      { 
        Jhat = sqrt(2*J+1); // the hat factor of J
      }
      t_start_omp = omp_get_wtime(); // profiling (s)
      #pragma omp parallel for schedule(dynamic,1) // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab|
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          int na = oa.n; // this is just...
          int nb = ob.n;
          int nc = oc.n;
          int nd = od.n;
          int la = oa.l;
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
          int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          for (int S=0; S<=1; S++) // sum over total spin...
          {
            int Seval;
            Seval = 2*S*(S + 1) - 3;
            int L_min = std::max(abs(la-lb),abs(lc-ld));
            int L_max = std::min(la+lb,lc+ld);
            for (int L = L_min; L <= L_max; L++) // ...and sum over orbital angular momentum
            {
              double sumMT = 0; // for the Moshinsky transformation
              double sumMTas = 0; // (anti-symmetric part)              
              double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
              double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
              double nNJab = normab*AngMom::NineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
              double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
              double nNJcd = normcd*AngMom::NineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
              double nNJdc = normcd*AngMom::NineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
              double bulk = Seval*nNJab*nNJcd; // bulk product of the above
              double bulkas = Seval*nNJab*nNJdc; // (anti-symmetric part)
              int tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
              for (int nr = 0; nr <= tempmaxnr; nr++)
              {
                double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of my thesis
                if ((npr >= 0) and ((eps_cd-eps_ab)%2 == 0))
                {
                  int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                  for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                  {
                    int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                    int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                    for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                    {
                      int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of my thesis
                      double integral = 0;
                      double normJrel;
                      double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                      double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                      double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
                      int minJrel= abs(lr-S);
                      int maxJrel = lr+S;
                      for (int Jrel = minJrel; Jrel<=maxJrel; Jrel++)
                      {
                        normJrel = sqrt((2*Jrel+1)*(2*L+1))*phase(L+lr+S+J)*AngMom::SixJ(Lam,lr,L,S,J,Jrel);
                        integral += normJrel*normJrel*GetM0nuIntegral(e2max,nr,lr,npr,lr,Jrel,hw,transition,Eclosure,src,IntList); // grab the pre-calculated integral wrt dq and dr from the IntList of the modelspace class
                      }//end of for-loop over Jrel
                      sumMT += Df*Di*integral; // perform the Moshinsky transformation
                      sumMTas += Df*asDi*integral; // (anti-symmetric part)
                    } // end of for-loop over: lr
                  } // end of for-loop over: Ncom
                } // end of if: npr \in \Nat_0
              } // end of for-loop over: nr
              sumLS += bulk*sumMT; // perform the LS-coupling sum
              sumLSas += bulkas*sumMTas; // (anti-symmetric part)
            } // end of for-loop over: L
          } // end of for-loop over: S
          double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          M0nuGT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
      M0nuGT_TBME.profiler.timer["M0nuGT_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
    } // end of for-loop over: auto
    std::cout<<"...done calculating M0nu TBMEs"<<std::endl;
    M0nuGT_TBME.profiler.timer["M0nuGT_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    M0nuGT_TBME.profiler.timer["M0nuGT_Op"] += omp_get_wtime() - t_start; // profiling (r)
    return M0nuGT_TBME;
  }