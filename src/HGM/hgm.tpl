//  ******************************************************************
//  Halibut Growth Model
//
//  Created by Martell on 2012-10-04.
//  Copyright (c) 2012. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION
	// Model dimensions
	init_int nbin;
	init_int nage;
	init_int narea;
	init_int nsex;
	
	init_vector ilen(1,nbin);
	init_vector jage(1,nage);
	int a1;
	int a2;
	LOC_CALCS
		a1 = min(jage);
		a2 = max(jage);
	END_CALCS
	
	
	init_4darray data(1,nsex,1,narea,1,nbin,1,nage);
	
	
	

PARAMETER_SECTION
	// Hyperparameters for growth model.
	init_vector   log_mu_l1(1,nsex,1);
	init_vector   log_mu_l2(1,nsex,1);
	init_vector  log_mu_rho(1,nsex,1);
	init_vector log_mu_lam1(1,nsex,1);
	init_vector log_mu_lam2(1,nsex,1);
	
	
	init_bounded_vector  log_mu_lhat(1,nsex,-5,5,2);
	init_bounded_vector  log_sd_lhat(1,nsex,-5,5,2);
	
	
	init_bounded_matrix   log_l1_dev(1,nsex,1,narea,-5.0,5.0,2);
	init_bounded_matrix   log_l2_dev(1,nsex,1,narea,-5.0,5.0,2);
	init_bounded_matrix  log_rho_dev(1,nsex,1,narea,-5.0,5.0,3);
	init_bounded_matrix log_lam1_dev(1,nsex,1,narea,-5.0,5.0,3);
	init_bounded_matrix log_lam2_dev(1,nsex,1,narea,-5.0,5.0,3);
	init_bounded_matrix log_lhat_dev(1,nsex,1,narea,-5.0,5.0,4);
	init_bounded_matrix log_sdlh_dev(1,nsex,1,narea,-5.0,5.0,4);
	init_bounded_matrix        log_Z(1,nsex,1,narea,-5.0,5.0,-3);
	
	
	objective_function_value f;
	
	number   fpen;
	
	matrix     l1(1,nsex,1,narea);
	matrix     l2(1,nsex,1,narea);
	matrix    rho(1,nsex,1,narea);
	matrix   lam1(1,nsex,1,narea);
	matrix   lam2(1,nsex,1,narea);
	matrix   lhat(1,nsex,1,narea);
	matrix   sdlh(1,nsex,1,narea);
	matrix      Z(1,nsex,1,narea);
	matrix   lvec(1,nsex,1,narea);
	matrix   pvec(1,nsex,1,narea);
	
	3darray mu_la(1,nsex,1,narea,1,nage);
	3darray sd_la(1,nsex,1,narea,1,nage);
	3darray    pa(1,nsex,1,narea,1,nage);   // age-proportions
	3darray    sl(1,nsex,1,narea,1,nbin);	// size-selectivity 
	
	
	4darray P(1,nsex,1,narea,1,nbin,1,nage);
	4darray V(1,nsex,1,narea,1,nbin,1,nage); // predicted proportions at length|age

INITIALIZATION_SECTION
	log_mu_l1    3.912023;
	log_mu_l2    5.010635;
	log_mu_rho  -0.20;
	log_mu_lhat  4.40;
	log_sd_lhat  2.25;
	
	log_mu_lam1  1.20;
	log_mu_lam2  0.40;
	      log_Z -1.61;

PROCEDURE_SECTION

	/*MAIN*/
	calcLengthAtAge();
	calcLengthAgeTransition();
	calcSizeSelectivity();
	calcAgeProportions();
	calcPredictedObservations();
	calcObjectiveFunction();
	cout<<f<<endl;
	/*    */

FUNCTION calcObjectiveFunction
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

FUNCTION calcPredictedObservations
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

FUNCTION calcSizeSelectivity
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

FUNCTION calcAgeProportions
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
	

FUNCTION calcLengthAtAge
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

FUNCTION calcLengthAgeTransition
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
	

FUNCTION dvar_vector dev_vector(dvar_vector &x, dvariable &pen)
  {
	/*A little tweek to make a dvar_vector sum to 0.*/
	dvariable s = mean(x);
	pen += 10000.0*s*s;
	x -=s;
	return(x);
  }

FUNCTION dvariable dnbinom(const dvector& x, const dvar_vector& mu, const double& k)
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

FUNCTION dvariable dnorm( const dvariable& x, const prevariable& mu, const prevariable& std )
  {

	if( std<=0 ) 
	{
		cerr<<"Standard deviation is less than or equal to zero in "
		"dnorm(const dvariable& x, const double& mu, const double& std)\n";
		return 0;
	}

	return 0.5*log(2.*M_PI)+log(std)+0.5*square(x-mu)/(std*std);
  }


FUNCTION dvariable dnorm( const dvariable& x, const prevariable& mu, const double& std )
  {

	if( std<=0 ) 
	{
		cerr<<"Standard deviation is less than or equal to zero in "
		"dnorm(const dvariable& x, const double& mu, const double& std)\n";
		return 0;
	}

	return 0.5*log(2.*M_PI)+log(std)+0.5*square(x-mu)/(std*std);
  }


REPORT_SECTION
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

TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
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
	
FINAL_SECTION
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

