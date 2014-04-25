	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#define TINY 1.e-12
	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <HGM.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nbin.allocate("nbin");
  nage.allocate("nage");
  narea.allocate("narea");
  nsex.allocate("nsex");
  ilen.allocate(1,nbin,"ilen");
  jage.allocate(1,nage,"jage");
		a1 = min(jage);
		a2 = max(jage);
  data.allocate(1,nsex,1,narea,1,nbin,1,nage,"data");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_mu_l1.allocate(1,nsex,1,"log_mu_l1");
  log_mu_l2.allocate(1,nsex,1,"log_mu_l2");
  log_mu_rho.allocate(1,nsex,1,"log_mu_rho");
  log_mu_lam1.allocate(1,nsex,1,"log_mu_lam1");
  log_mu_lam2.allocate(1,nsex,1,"log_mu_lam2");
  log_mu_lhat.allocate(1,nsex,-5,5,2,"log_mu_lhat");
  log_sd_lhat.allocate(1,nsex,-5,5,2,"log_sd_lhat");
  log_l1_dev.allocate(1,nsex,1,narea,-5.0,5.0,2,"log_l1_dev");
  log_l2_dev.allocate(1,nsex,1,narea,-5.0,5.0,2,"log_l2_dev");
  log_rho_dev.allocate(1,nsex,1,narea,-5.0,5.0,3,"log_rho_dev");
  log_lam1_dev.allocate(1,nsex,1,narea,-5.0,5.0,3,"log_lam1_dev");
  log_lam2_dev.allocate(1,nsex,1,narea,-5.0,5.0,3,"log_lam2_dev");
  log_lhat_dev.allocate(1,nsex,1,narea,-5.0,5.0,4,"log_lhat_dev");
  log_sdlh_dev.allocate(1,nsex,1,narea,-5.0,5.0,4,"log_sdlh_dev");
  log_Z.allocate(1,nsex,1,narea,-5.0,5.0,-3,"log_Z");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  l1.allocate(1,nsex,1,narea,"l1");
  #ifndef NO_AD_INITIALIZE
    l1.initialize();
  #endif
  l2.allocate(1,nsex,1,narea,"l2");
  #ifndef NO_AD_INITIALIZE
    l2.initialize();
  #endif
  rho.allocate(1,nsex,1,narea,"rho");
  #ifndef NO_AD_INITIALIZE
    rho.initialize();
  #endif
  lam1.allocate(1,nsex,1,narea,"lam1");
  #ifndef NO_AD_INITIALIZE
    lam1.initialize();
  #endif
  lam2.allocate(1,nsex,1,narea,"lam2");
  #ifndef NO_AD_INITIALIZE
    lam2.initialize();
  #endif
  lhat.allocate(1,nsex,1,narea,"lhat");
  #ifndef NO_AD_INITIALIZE
    lhat.initialize();
  #endif
  sdlh.allocate(1,nsex,1,narea,"sdlh");
  #ifndef NO_AD_INITIALIZE
    sdlh.initialize();
  #endif
  Z.allocate(1,nsex,1,narea,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  lvec.allocate(1,nsex,1,narea,"lvec");
  #ifndef NO_AD_INITIALIZE
    lvec.initialize();
  #endif
  pvec.allocate(1,nsex,1,narea,"pvec");
  #ifndef NO_AD_INITIALIZE
    pvec.initialize();
  #endif
  mu_la.allocate(1,nsex,1,narea,1,nage,"mu_la");
  #ifndef NO_AD_INITIALIZE
    mu_la.initialize();
  #endif
  sd_la.allocate(1,nsex,1,narea,1,nage,"sd_la");
  #ifndef NO_AD_INITIALIZE
    sd_la.initialize();
  #endif
  pa.allocate(1,nsex,1,narea,1,nage,"pa");
  #ifndef NO_AD_INITIALIZE
    pa.initialize();
  #endif
  sl.allocate(1,nsex,1,narea,1,nbin,"sl");
  #ifndef NO_AD_INITIALIZE
    sl.initialize();
  #endif
  P.allocate(1,nsex,1,narea,1,nbin,1,nage,"P");
  #ifndef NO_AD_INITIALIZE
    P.initialize();
  #endif
  V.allocate(1,nsex,1,narea,1,nbin,1,nage,"V");
  #ifndef NO_AD_INITIALIZE
    V.initialize();
  #endif
}

void model_parameters::initializationfunction(void)
{
  log_mu_l1.set_initial_value(3.912023);
  log_mu_l2.set_initial_value(5.010635);
  log_mu_rho.set_initial_value(-0.20);
  log_mu_lhat.set_initial_value(4.40);
  log_sd_lhat.set_initial_value(2.25);
  log_mu_lam1.set_initial_value(1.20);
  log_mu_lam2.set_initial_value(0.40);
  log_Z.set_initial_value(-1.61);
}

void model_parameters::userfunction(void)
{
  f =0.0;
	/*MAIN*/
	calcLengthAtAge();
	calcLengthAgeTransition();
	calcSizeSelectivity();
	calcAgeProportions();
	calcPredictedObservations();
	calcObjectiveFunction();
	cout<<f<<endl;
	/*    */
}

void model_parameters::calcObjectiveFunction(void)
{
	/**
		Compute the objective function (f) for minimization
		Using the negative log of the multinomial distribution as the likelihood.
	*/
	lvec.initialize();
	int i,j,k,l;
	double tau2 = 8.0;
	for(i=1;i<=nsex;i++)
	{
		for(j=1;j<=narea;j++)
		{
			for(l=1;l<=nbin;l++)
			{
				lvec(i)(j) += dmultinom( data(i)(j)(l), V(i)(j)(l)+TINY );
				//lvec(i)(j) += dnbinom( data(i)(j)(l), V(i)(j)(l), 8.0 );
			}
		}
	}
	/*priors*/
	pvec.initialize();
	for(i=1;i<=nsex;i++)
	{
		for(j=1;j<=narea;j++)
		{
			pvec(i)(j)  = dnorm(log(l1(i)(j)),log_mu_l1(i),0.2);
			pvec(i)(j) += dnorm(log(l2(i)(j)),log_mu_l2(i),0.2);
			pvec(i)(j) += dnorm(log(rho(i)(j)),log_mu_rho(i),0.2);
		}
	}
	f = sum(lvec) + sum(pvec) + fpen;
}

void model_parameters::calcPredictedObservations(void)
{
	/*
		Compute the predicted length-age matrixes
		V(i,j,k,k) = pa(i,j,k)*sa(i,j,l)*P(i,j,l,k)
	*/
	V.initialize();
	int i,j,k,l;
	for(i=1;i<=nsex;i++)
	{
		for(j=1;j<=narea;j++)
		{
			for(l=1;l<=nbin;l++)
			{
				V(i)(j)(l) = elem_prod(pa(i)(j)*sl(i)(j)(l),P(i)(j)(l));
			}
			V(i)(j) /= sum(V(i)(j));
		}
	}
}

void model_parameters::calcSizeSelectivity(void)
{
	/*
		Compute size-based selectivity using a 
		simple logistic function.
		Indices:
			i = sex
			j = area
	*/
	sl.initialize();
	int i,j;
	double bw2 = 0.5*(ilen(2)-ilen(1));
	for(i=1;i<=nsex;i++)
	{
		lhat(i) = mfexp(log_mu_lhat(i) + dev_vector(log_lhat_dev(i),fpen));
		sdlh(i) = mfexp(log_sd_lhat(i) + dev_vector(log_sdlh_dev(i),fpen))+TINY;
		for(j=1;j<=narea;j++)
		{
			sl(i)(j) = plogis(ilen+bw2,lhat(i)(j),sdlh(i)(j));
		}
	}
}

void model_parameters::calcAgeProportions(void)
{
	/*
		Compute the relative proporitons-at-age based on
		estimates of Z in each of the areas
		TRICKY PART:  
		Need to compute selectivity-at-age which
		is a function of size-at-age, variation in size-
		at-age and size-selectivity.
		Indices:
			i = sex
			j = area
			k = age
	*/
	pa.initialize();
	int i,j,k;
	dvar_vector sa(1,nage);
	for(i=1;i<=nsex;i++)
	{
		Z(i) = mfexp(log_Z(i));
		for(j=1;j<=narea;j++)
		{
			/*Age-specific selectivity (intgrated over size-at-age)*/
			sa = trans(P(i)(j))*sl(i)(j);
			pa(i)(j)(1) = 1.0;
			for(k=2;k<=nage;k++)
			{
				pa(i)(j)(k) = pa(i)(j)(k-1) * mfexp(-Z(i)(j)*sa(k-1));
				if(k==nage)
				{
					pa(i)(j)(k) /= (1.0 - mfexp(-Z(i)(j)*sa(k)));
				}
			}
		}
	}
}

void model_parameters::calcLengthAtAge(void)
{
	/*
		Indicies:
			i = sex
			j = area
			k = age
	*/
	int i,j,k;
	dvar_vector t1(1,nage);
	for(i=1;i<=nsex;i++)
	{
		l1(i)   = mfexp(log_mu_l1(i)   +   dev_vector(log_l1_dev(i),fpen));
		l2(i)   = mfexp(log_mu_l2(i)   +   dev_vector(log_l2_dev(i),fpen));
		rho(i)  = mfexp(log_mu_rho(i)  +  dev_vector(log_rho_dev(i),fpen));
		lam1(i) = mfexp(log_mu_lam1(i) + dev_vector(log_lam1_dev(i),fpen));
		lam2(i) = mfexp(log_mu_lam2(i) + dev_vector(log_lam2_dev(i),fpen));
		for(j=1;j<=narea;j++)
		{
			t1          = 1.0-pow(rho(i)(j),jage-double(a1));
			t1         /= 1.0-pow(rho(i)(j),double(a2)-double(a1));
			mu_la(i)(j) = l1(i)(j)+(l2(i)(j)-l1(i)(j))*t1;
			sd_la(i)(j) = lam1(i)(j)*exp( lam2(i)(j)*(-1.0+2.0*t1) );
		}
	}
}

void model_parameters::calcLengthAgeTransition(void)
{
	/*
		This function calculates the P(l|a) matrix, which is read as
		the probability of being in length interval (l) given age (a).
		Indices:
			i = sex
			j = area
			k = age
			l = length
	*/
	int i,j,k,l;
	double bw = ilen(2)-ilen(1);
	dvariable z1;
	dvariable z2;
	dvar_matrix ptrans(1,nage,1,nbin);
	for(i=1;i<=nsex;i++)
	{
		for(j=1;j<=narea;j++)
		{
			for(k=1;k<=nage;k++)
			{
				for(l=1;l<=nbin;l++)
				{
					z1 = ( ilen(l)    - mu_la(i)(j)(k) ) / sd_la(i)(j)(k);
					z2 = ( ilen(l)+bw - mu_la(i)(j)(k) ) / sd_la(i)(j)(k);
					ptrans(k,l) = cumd_norm(z2) - cumd_norm(z1);
				}
				ptrans(k) /= sum(ptrans(k));
			}
			P(i)(j) = trans(ptrans);
		}
	}
	//ofstream log("hgmP.log");
	//log<<P<<endl;
	//exit(1);
}

dvar_vector model_parameters::dev_vector(dvar_vector &x, dvariable &pen)
{
  {
	/*A little tweek to make a dvar_vector sum to 0.*/
	dvariable s = mean(x);
	pen += 10000.0*s*s;
	x -=s;
	return(x);
  }
}

dvariable model_parameters::dnbinom(const dvector& x, const dvar_vector& mu, const double& k)
{
  {
	//the observed counts are in x
	//mu is the predicted mean
	//k is the overdispersion parameter
	if (k<0.0)
	{
		cerr<<"k is <=0.0 in dnbinom()";
		return(0.0);
	}
	RETURN_ARRAYS_INCREMENT();
	int i,imin,imax;
	imin=x.indexmin();
	imax=x.indexmax();
	dvariable loglike = 0.;
	for(i = imin; i<=imax; i++)
	{
		//cout<<gammln(k+x(i))-gammln(k)-gammln(x(i)+1.)+k*log(k)-k*log(mu(i)+k)+x(i)*log(mu(i)+TINY)-x(i)*log(mu(i)+k)<<endl;
		loglike += gammln(k+x(i))-gammln(k)-gammln(x(i)+1.)+k*log(k)-k*log(mu(i)+k)+x(i)*log(mu(i)+TINY)-x(i)*log(mu(i)+k);
	}
	RETURN_ARRAYS_DECREMENT();
	return(-loglike);
  }
}

dvariable model_parameters::dnorm( const dvariable& x, const prevariable& mu, const prevariable& std )
{
  {
	if( std<=0 ) 
	{
		cerr<<"Standard deviation is less than or equal to zero in "
		"dnorm(const dvariable& x, const double& mu, const double& std)\n";
		return 0;
	}
	return 0.5*log(2.*M_PI)+log(std)+0.5*square(x-mu)/(std*std);
  }
}

dvariable model_parameters::dnorm( const dvariable& x, const prevariable& mu, const double& std )
{
  {
	if( std<=0 ) 
	{
		cerr<<"Standard deviation is less than or equal to zero in "
		"dnorm(const dvariable& x, const double& mu, const double& std)\n";
		return 0;
	}
	return 0.5*log(2.*M_PI)+log(std)+0.5*square(x-mu)/(std*std);
  }
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
	REPORT(nbin  );
	REPORT(nage  );
	REPORT(narea );
	REPORT(nsex  );
	REPORT(jage  ); 
	REPORT(ilen  );
	REPORT(mu_la);
	REPORT(sd_la);
	REPORT(sl);
	report<<"P\n"<<P(1)(1)<<endl;
	report<<"V\n"<<V(1)(1)<<endl;
	REPORT(data);
}

void model_parameters::final_calcs()
{
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
