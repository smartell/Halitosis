//  ******************************************************************
//  vonBH
//  
//  Created by Martell on 2012-09-28.
//  Copyright (c) 2012. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION
	init_int n;
	init_int narea;
	init_int nsex;
	init_matrix data(1,n,1,4);
	
	//#area	 age	 fl	 sex
	ivector area(1,n);
	vector   age(1,n);
	vector    fl(1,n);
	ivector  sex(1,n);

	LOC_CALCS
		area = ivector(column(data,1));
		age  = column(data,2);
		fl   = column(data,3);
		sex  = ivector(column(data,4));
	END_CALCS
	
	// referenece points
	matrix f01(1,nsex,1,narea);

PARAMETER_SECTION
	init_vector log_mu_linf(1,nsex);
	init_vector  log_mu_vbk(1,nsex);
	init_vector       mu_to(1,nsex,-1);
	init_vector  log_mu_sig(1,nsex);
	
	init_vector log_sd_linf(1,nsex,-4);
	init_vector  log_sd_vbk(1,nsex,-4);
	init_vector   log_sd_to(1,nsex,-4);
	init_vector  log_sd_sig(1,nsex,-4);
	
	init_bounded_matrix linf_dev(1,nsex,1,narea,-5,5,2);
	init_bounded_matrix    k_dev(1,nsex,1,narea,-5,5,2);
	init_bounded_matrix   to_dev(1,nsex,1,narea,-5,5,2);
	init_bounded_matrix  sig_dev(1,nsex,1,narea,-5,5,3);
	
	objective_function_value f;
	
	number fpen;
	
	matrix linf(1,nsex,1,narea);
	matrix  vbk(1,nsex,1,narea);
	matrix   to(1,nsex,1,narea);
	matrix  sig(1,nsex,1,narea);
	matrix lvec(1,nsex,1,narea);
	matrix pvec(1,nsex,1,narea);

	vector fl_hat(1,n);
	vector epsilon(1,n);
	sdreport_number sd_linf;
	
INITIALIZATION_SECTION
	log_mu_linf 5.0;
	log_mu_vbk -1.7;
	    mu_to  -0.5;
	log_mu_sig -1.7;
	
	log_sd_linf 0.0;
	log_sd_vbk -1.0;
	log_sd_to  -1.6;
	log_sd_sig  0.0;
	
	
PROCEDURE_SECTION
	vonb_model();
	sd_linf = mfexp(log_mu_linf(1));
	calc_objective_function();
	if(mceval_phase())
		mcmcReport();
	
FUNCTION vonb_model
	// i = observation
	// j = sex
	// k = area
	int i,j,k;
	
	/*Initialize parameters*/
	fpen = 0;
	for(j=1;j<=nsex;j++)
	{
		linf(j) = mfexp(log_mu_linf(j) + dev_vector(linf_dev(j),fpen));
		 vbk(j) = mfexp(log_mu_vbk(j) + dev_vector(k_dev(j),fpen));
		  to(j) = (mu_to(j) + dev_vector(to_dev(j),fpen));
		 sig(j) = mfexp(log_mu_sig(j) + dev_vector(sig_dev(j),fpen));
	}
	
	/*Compute predicted lengths*/
	dvariable tmpl;
	for(i=1;i<=n;i++)
	{
		j    = sex(i);
		k    = area(i);
		tmpl = linf(j,k) * (1.0 - mfexp(-vbk(j,k)*(age(i)-to(j,k))));
		fl_hat(i) = tmpl;
		epsilon(i) = (fl(i) - fl_hat(i))/sig(j,k);
	}

FUNCTION calc_objective_function
	lvec.initialize();
	pvec.initialize();
	int i,j,k;
	
	/*negative log likelihoods*/
	for(i=1;i<=n;i++)
	{
		j         = sex(i);
		k         = area(i);
		lvec(j,k)+= dnorm(fl(i),fl_hat(i),sig(j,k));
	}
	
	/*prior densities*/
	for(j=1;j<=nsex;j++)
	{
		for(k=1;k<=narea;k++)
		{
			pvec(j,k)  = dnorm(log(linf(j,k)),log_mu_linf(j),mfexp(log_sd_linf(j)));
			pvec(j,k) += dnorm(log(vbk(j,k)),log_mu_vbk(j),mfexp(log_sd_vbk(j)));
			pvec(j,k) += dnorm((to(j,k)),mu_to(j),mfexp(log_sd_to(j)));
			pvec(j,k) += dnorm(log(sig(j,k)),log_mu_sig(j),mfexp(log_sd_sig(j)));
			
		}
	}
	
	f = sum(lvec) + sum(pvec) + fpen;
	
	
FUNCTION dvar_vector dev_vector(dvar_vector &x, dvariable &pen)
  {
	/*A little tweek to make a dvar_vector sum to 0.*/
	dvariable s = mean(x);
	pen += 10000.0*s*s;
	x -=s;
	return(x);
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

REPORT_SECTION
	REPORT(area);
	REPORT(sex);
	REPORT(age);
	REPORT(fl);
	REPORT(fl_hat);
	REPORT(epsilon);
	
	int max_age = max(age);
	dvector iage(1,max_age);
	iage.fill_seqadd(1,1);
	dmatrix p_fl(1,nsex,1,max_age);
	for(int j=1;j<=nsex;j++)
	{
		double linf = value(exp( log_mu_linf(j) ));
		double vbk  = value(exp( log_mu_vbk(j)  ));
		double to   = value(mu_to(j));
		p_fl(j) = linf*(1.-exp(-vbk*(iage-to)));
	}
	REPORT(iage);
	REPORT(p_fl);
	
	//dvector linf(1,nsex);
	//dvector  vbk(1,nsex);
	//dvector   to(1,nsex);
	//linf = value(exp(log_mu_linf));
	// vbk = value(exp(log_mu_vbk));
	//  to = value(mu_to);
	
	REPORT(linf);
	REPORT(vbk);
	REPORT(to);
	
	
FUNCTION mcmcReport
  {
	/*Write MCMC report*/
	static int nf=0;
	
	if(nf==0)
	{
		cout<<"Writing MCMC report... Please wait"<<endl;
		ofstream ofs("vonBH.mcmc");
		ofs<<"sex \t area \t F0.1"<<endl;
	}
	nf ++;
	cout<<nf<<endl;
	ofstream ofs("vonBH.mcmc",ios::app);
	
	calcReferencePoints();
	
	int j,k;
	for(j=1;j<=nsex;j++)
	{
		for(k=1;k<=narea;k++)
		{
			ofs<<j<<"\t"<<k<<"\t"<<f01(j,k)<<endl;
		}
		
	}
	
  }
	
FUNCTION calcReferencePoints
  {
	/*
	Calculate Yield per Recruit based reference points based on growth
	in each regulatory area, assuming coast-wide length-based selectivity
	from the 2012 assessment model.
	
	Dependencies (ypr.cpp, and selex.cpp)
	*/
	int i,j,k;
	//double f01;
	int nage = max(age);
	dvector age(1,nage);
	age.fill_seqadd(1,1);
	
	/*Selectivities*/
	dvector bin(1,8);
	dmatrix CSelL(1,2,1,8);
	bin.fill("{50,60,70,80,90,100,110,120}");
	/*Based on results from Hare's 2012 assessment*/
	CSelL(1).fill("{1.63180e-09,3.25740e-09,0.0603022,0.300891,0.630344,0.913893,1,1}");
	CSelL(2).fill("{2.17159e-09,3.96878e-03,0.0567109,0.281436,0.585461,0.835614,1,1}");
	Selex c_selexfun;
	
	/*loop over sexs and areas and calculate F0.1*/
	double a = 6.82e-6;
	double b = 3.24;
	dvector m(1,2);
	m(1) = 0.15;	//female m
	m(2) = 0.135;	//male m
	
	
	dvector la(1,nage);
	dvector wa(1,nage);
	dvector sa(1,nage);
	for(j=1;j<=nsex;j++)
	{
		for(k=1;k<=narea;k++)
		{
			/*Growth*/
			double m_linf = value(linf(j,k));
			double m_vbk  = value(vbk(j,k) );
			double m_to   = value(to(j,k)  );
			la = m_linf*(1.0-exp(-m_vbk*(age-m_to)));
			wa = a * pow(la,b);
			
			/*Selectivity*/
			sa = c_selexfun.linapprox(bin,CSelL(j),la);
			
			/*F0.1 reference point*/
			yieldPerRecruit c_ypr(m(j),wa,sa);
			f01(j,k) = c_ypr.getf01();
			if(!mceval_phase())cout<<"Sex ="<<j<<" Area = "<<k<<" F0.1="<<f01(j,k)<<endl;
		}
	}
	
  }
	
FUNCTION double ypr(const double& fe, const int& id)
  {
	/*Calculate the yield per recruit for a given fishing rate fe*/
	/*
		Algorithm
		ypr = fe * sum_i(sum_j(lz_{ij}*wa_{ij}*va_{ij})) 
		
		index: i= sex, j= age
		1) Compute sex-specific length-at-age
		2) Compute sex-specific weight-at-age
		3) Compute sex-specific selectivity-at-age
		4) Compute sex-specific survivorship
	*/
	int i,j;
	int nage = max(age);
	dvector age(1,nage);
	age.fill_seqadd(1,1);
	dvector m(1,2);
	m(1) = 0.15;	//female m
	m(2) = 0.135;	//male m
	double a = 6.82e-6;
	double b = 3.24;
	dmatrix la(1,2,1,nage);
	dmatrix wa(1,2,1,nage);
	dmatrix lx(1,2,1,nage);
	dmatrix lz(1,2,1,nage);
	dmatrix sx(1,2,1,nage);
	Selex cselex;
	dvector bin(1,8);
	bin.fill("{50,60,70,80,90,100,110,120}");
	dmatrix CSelL(1,2,1,8);
	CSelL(1).fill("{1.63180e-09,3.25740e-09,0.0603022,0.300891,0.630344,0.913893,1,1}");
	CSelL(2).fill("{2.17159e-09,3.96878e-03,0.0567109,0.281436,0.585461,0.835614,1,1}");
	
	/*Growth*/
	la.initialize();
	wa.initialize();
	for(i=1;i<=2;i++)
	{
		la(i) = value(linf(i,id)*(1.0-exp(-vbk(i,id)*(age-to(i,id)))));
		wa(i) = a * pow(la(i),b);
		sx(i) = cselex.linapprox(bin,CSelL(i),la(i));
	}
	
	//	yieldPerRecruit(double fe, double m, dvector wa, dvector sa);
	yieldPerRecruit c_ypr(m(1),wa(1),sx(1));
	//cout<<"YPR="<<c_ypr.getypr(0.01)<<endl;
	//cout<<"YPR="<<c_ypr.getypr(0.2)<<endl;
	//cout<<"YPR="<<c_ypr.getypr(0.3)<<endl;
	
	
	/*Survivorship*/
	lx.initialize();
	lz.initialize();
	lx = 1.0;
	lz = 1.0;
	for(i=1;i<=2;i++)
	{
		for(j=2;j<=nage;j++)
		{
			lx(i,j) = lx(i,j-1)*exp(-m(i));
			lz(i,j) = lz(i,j-1)*exp(-m(i)-fe*sx(i,j));
		}
	}
	double ypr0 = c_ypr.getypr(0);
	//double ypr  = fe*sum(elem_prod(lz,elem_prod(wa,sx)));
	double ypr = c_ypr.getypr(fe);
	double f01 = c_ypr.getf01();
	cout<<fe<<"\t"<<ypr<<"\t"<<ypr0<<"\t"<<f01<<endl;
	return(ypr);
  }

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

	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	#include <selex.cpp>
	#include <ypr.cpp>
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
	calcReferencePoints();
	

