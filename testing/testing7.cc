Operator Tensor(ModelSpace& modelspace, double Eclosure, std::string src)
{
  bool reduced = true;
  double t_start, t_start_tbme, t_start_omp; // profiling (v)
  t_start = omp_get_wtime(); // profiling (s)
  std::string transition = "T";
  // run through the initial set-up routine
  double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
  int e2max = modelspace.GetE2max(); // 2*emax
  Operator M0nuT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
  std::cout<<"     reduced            =  "<<reduced<<std::endl;
  M0nuT_TBME.SetHermitian(); // it should be Hermitian
  int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
  const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
  const double prefact = Rnuc/(PI*PI); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
  modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
  std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency
  M0nuT_TBME.profiler.timer["M0nuT_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
  // create the TBMEs of M0nu
  // auto loops over the TBME channels and such
  std::cout<<"calculating M0nu TBMEs..."<<std::endl;
  t_start_tbme = omp_get_wtime(); // profiling (s)
  for (auto& itmat : M0nuT_TBME.TwoBody.MatEl)
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
    else //if (reduced == "R")
    {
      Jhat = sqrt(2*J + 1); // the hat factor of J
    }
    t_start_omp = omp_get_wtime(); // profiling (s)
    #pragma omp parallel for schedule(dynamic,1) // need to do: PreCalculateMoshinsky(), PreCalcT6j, and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
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
        double ja = oa.j2/2.0;
        double jb = ob.j2/2.0;
        double jc = oc.j2/2.0;
        double jd = od.j2/2.0; // ...for convenience
        int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
        int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
        double sumLS = 0; // for the wave functions decomposition
        double sumLSas = 0; // (anti-symmetric part)
        int S = 1;
        int Lf_min = std::max(std::abs(la-lb), std::abs(J-S));
        int Lf_max = std::min(la+lb, J+S);
        for (int Lf = Lf_min; Lf<=Lf_max; Lf++) // sum over angular momentum coupled to l_a and l_b
        {
          double normab = sqrt((2*Lf+1)*(2*S+1)*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
          double nNJab = normab*AngMom::NineJ(la,lb,Lf,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
          int Li_min = std::max(std::abs(lc-ld), std::max(std::abs(J-S), std::abs(Lf-2)));
          int Li_max = std::min( lc+ld, std::min( J+S, Lf+2) );
          for (int Li = Li_min; Li <= Li_max; Li++) // sum over angular momentum coupled to l_c and l_d
          { 
            double sumMT = 0; // for the Moshinsky transformation
            double sumMTas = 0; // (anti-symmetric part)
            double normcd = sqrt((2*Li+1)*(2*S+1)*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
            double nNJcd = normcd*AngMom::NineJ(lc,ld,Li,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
            double nNJdc = normcd*AngMom::NineJ(ld,lc,Li,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)

            double bulk = nNJab*nNJcd; // bulk product of the above
            double bulkas = nNJab*nNJdc; // (anti-symmetric part)
            for ( int lr=1; lr<=eps_ab; lr++)
            {
              for (int nr=0; nr<=(eps_ab-lr)/2; nr++)
              {
                int tempmaxNcom = std::min((eps_ab-2*nr-lr)/2, eps_cd);
                for (int Ncom=0; Ncom<=tempmaxNcom; Ncom++ )
                {
                  int Lam = eps_ab - 2*nr - lr - 2*Ncom;
                  if ( (Lam+2*Ncom) > eps_cd ) continue;
                  if ( (std::abs(Lam-lr)>Lf)  or ( (Lam+lr)<Lf) ) continue;
                  if ( (lr+Lam+eps_ab)%2>0 ) continue;
                  for (int npr=0; npr<=(eps_cd-2*Ncom-Lam)/2; npr++)
                  {
                      int lpr = eps_cd-2*Ncom-Lam-2*npr;
                      if (  (lpr+lr)%2 >0 ) continue;
                      if (lpr<1) continue;
                      if ( (std::abs(Lam-lpr)>Li)  or ( (Lam+lpr)<Li) ) continue;
                      double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,Lf); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                      double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nc,lc,nd,ld,Li); // " " " "
                      double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nd,ld,nc,lc,Li);// (anti-symmetric part)
                      double integral = 0;
                      double normJrel, normJrelp;
                      int minJrel = std::max(abs(lr-S),abs(lpr-S));
                      int maxJrel = std::min(lr+S,lpr+S);
                      for (int Jrel = minJrel; Jrel<=maxJrel; Jrel++)
                      {
                        if ( (std::abs(J-Jrel)>Lam)  or ( (Jrel+J)<Lam) ) continue;
                        normJrel  = sqrt((2*Jrel+1)*(2*Lf+1))*phase(Lf+lr+J+S)*AngMom::SixJ(Lam,lr,Lf,S,J,Jrel);
                        normJrelp = sqrt((2*Jrel+1)*(2*Li+1))*phase(Li+lpr+J+S)*AngMom::SixJ(Lam,lpr,Li,S,J,Jrel);
                        integral += normJrel*normJrelp*GetM0nuIntegral(e2max,nr,lr,npr,lpr,Jrel,hw,transition,Eclosure,src,IntList);
                      }
                      sumMT += Df*Di*integral; // perform the Moshinsky transformation
                      sumMTas += Df*asDi*integral; // (anti-symmetric part)
                    } // end of for-loop over: lpr
                  } // end of for-loop over: Ncom
                } // end of for-loop over: nr
              } // end of for-loop over: lr
            sumLS += bulk*sumMT; // perform the LS-coupling sum
            sumLSas += bulkas*sumMTas; // (anti-symmetric part)
          } // end of for-loop over: Li
        } // end of for-loop over: Lf        
        // double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*sumLS; // compute the final matrix element, anti-symmetrize          
        double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
        M0nuT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
      } // end of for-loop over: iket
    } // end of for-loop over: ibra
    M0nuT_TBME.profiler.timer["M0nuT_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
  } // end of for-loop over: auto
  M0nuT_TBME.profiler.timer["M0nuT_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
  M0nuT_TBME.profiler.timer["M0nuT _Op"] += omp_get_wtime() - t_start; // profiling (r)
  return M0nuT_TBME;
}