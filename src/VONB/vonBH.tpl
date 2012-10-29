//  ******************************************************************
//  vonBH
//  
//  Created by Martell on 2012-09-28.
//  Copyright (c) 2012. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION
	// referenece ages
	int age_1;
	int age_2;
	!! age_1 = 2;
	!! age_2 = 20;
	
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
	init_bounded_number   m(0,0.5);
	init_vector   log_mu_l1(1,nsex);
	init_vector   log_mu_l2(1,nsex);
	init_vector  log_mu_rho(1,nsex);
	init_vector    log_mu_b(1,nsex,2);
	init_vector   log_mu_cv(1,nsex,3);
	
	init_vector  log_sig_l1(1,nsex,-2);
	init_vector  log_sig_l2(1,nsex,-2);
	init_vector log_sig_rho(1,nsex,-2);
	init_vector   log_sig_b(1,nsex,-2);
	init_vector  log_sig_cv(1,nsex,-2);
	
	
	init_bounded_matrix  l1_dev(1,nsex,1,narea,-5,5,3);
	init_bounded_matrix  l2_dev(1,nsex,1,narea,-5,5,3);
	init_bounded_matrix rho_dev(1,nsex,1,narea,-5,5,3);
	init_bounded_matrix   b_dev(1,nsex,1,narea,-5,5,4);
	init_bounded_matrix  cv_dev(1,nsex,1,narea,-5,5,5);
	
	objective_function_value f;
	
	number fpen;
	
	matrix   l1(1,nsex,1,narea);
	matrix   l2(1,nsex,1,narea);
	matrix  rho(1,nsex,1,narea);
	matrix    b(1,nsex,1,narea);
	matrix   cv(1,nsex,1,narea);
	matrix lvec(1,nsex,1,narea);
	matrix pvec(1,nsex,1,narea);

	vector fl_hat(1,n);
	vector sd_fl(1,n);
	vector epsilon(1,n);
	sdreport_number sd_linf;
	
INITIALIZATION_SECTION
       log_mu_l1  2.00;
       log_mu_l2  4.90;
      log_mu_rho -0.15;
        log_mu_b  0.00;
       log_mu_cv -2.30;
    
      log_sig_l1  -15.24;
      log_sig_l2  -0.07;
     log_sig_rho  -5.30;
	   log_sig_b  -12.30;
	  log_sig_cv  -2.30;
	
PROCEDURE_SECTION
	vonb_model();
	sd_linf = mfexp(log_mu_l2(1));
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
		 l1(j) = mfexp( log_mu_l1(j) + dev_vector( l1_dev(j),fpen));
		 l2(j) = mfexp( log_mu_l2(j) + dev_vector( l2_dev(j),fpen));
		rho(j) = mfexp(log_mu_rho(j) + dev_vector(rho_dev(j),fpen));
		  b(j) = mfexp(  log_mu_b(j) + dev_vector(  b_dev(j),fpen));
		 cv(j) = mfexp( log_mu_cv(j) + dev_vector( cv_dev(j),fpen));
	}
	
	/*Compute predicted lengths fl_hat*/
	for(i=1;i<=n;i++)
	{
		j    = sex(i);
		k    = area(i);
		
		fl_hat(i) = get_len(age(i),age_1,age_2,l1(j,k),l2(j,k),rho(j,k),b(j,k));
		sd_fl(i)  =  cv(j,k)*fl_hat(i);
		epsilon(i) = (fl(i) - fl_hat(i))/sd_fl(i);
	}
	
FUNCTION dvariable get_len(const double& age,const int& age_1, const int& age_2, const dvariable& l1,const dvariable& l2,const dvariable& rho, const dvariable& b)
  {
		dvariable t1 = b*log(l1);
		dvariable t3 = mfexp(t1);
		dvariable t5 = mfexp(b*log(l2)) - t3;
		dvariable t6 = 1.0-mfexp((age-age_1)*log(rho));
		dvariable t7 = 1.0-mfexp((age_2-age_1)*log(rho));
		dvariable t9 = t6/t7;
		dvariable t0 = 1.0/b;
		
		dvariable len = mfexp(t0*(log(t3+t5*t9)));
		
		//dvariable t1 = pow(l1,b);
		//dvariable t3 = pow(l2,b);
		//dvariable t5 = 1.0-pow(rho,age-double(age_1));
		//dvariable t7 = 1.0-pow(rho,double(age_2)-double(age_1));
		//dvariable t9 = 1.0/b;
		//
		//dvariable len = pow(t1 + (t3-t1)*t5/t7,t9);
		return(len);
  }

FUNCTION double get_len(const double& age,const int& age_1, const int& age_2, const double& l1,const double& l2,const double& rho, const double& b)
  {
		double t1 = b*log(l1);
		double t3 = mfexp(t1);
		double t5 = mfexp(b*log(l2)) - t3;
		double t6 = 1.0-mfexp((age-age_1)*log(rho));
		double t7 = 1.0-mfexp((age_2-age_1)*log(rho));
		double t9 = t6/t7;
		double t0 = 1.0/b;
        
		double len = mfexp(t0*(log(t3+t5*t9)));

		//dvariable t1 = pow(l1,b);
		//dvariable t3 = pow(l2,b);
		//dvariable t5 = 1.0-pow(rho,age-double(age_1));
		//dvariable t7 = 1.0-pow(rho,double(age_2)-double(age_1));
		//dvariable t9 = 1.0/b;
		//
		//dvariable len = pow(t1 + (t3-t1)*t5/t7,t9);
		return(len);
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
		lvec(j,k)+= dnorm(fl(i),fl_hat(i),sd_fl(i));
	}
	
	/*prior densities*/
	for(j=1;j<=nsex;j++)
	{
		for(k=1;k<=narea;k++)
		{
			pvec(j,k)  = dnorm(log(l1(j,k)), log_mu_l1(j), mfexp(log_sig_l1(j)) );
			pvec(j,k) += dnorm(log(l2(j,k)), log_mu_l2(j), mfexp(log_sig_l2(j)) );
			pvec(j,k) += dnorm(log(rho(j,k)),log_mu_rho(j),mfexp(log_sig_rho(j)));
			pvec(j,k) += dnorm(log(b(j,k)),  log_mu_b(j),  mfexp(log_sig_b(j))  );
			pvec(j,k) += dnorm(log(cv(j,k)), log_mu_cv(j), mfexp(log_sig_cv(j)) );	
		}
	}
	
	/* prior for natural mortality rate */
	double a = 25;
	double s = 1./0.006;
	dvariable m_prior;
	m_prior  = dgamma(m,a,s);
	m_prior += dnorm(-m/log_mu_rho(1),double(1.5),double(0.1));
	
	f = sum(lvec) + sum(pvec) + fpen + m_prior;
	
	
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
FUNCTION dvariable dnorm( const dvariable& x, const double& mu, const double& std )
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
	
	dvector l1 = value(mfexp(log_mu_l1));
	dvector l2 = value(mfexp(log_mu_l2));
	dvector rho = value(mfexp(log_mu_rho));
	dvector b = value(mfexp(log_mu_b));
	dvector cv = value(mfexp(log_mu_cv));
	REPORT(l1);
	REPORT(l2);
	REPORT(rho);
	REPORT(b);
	REPORT(cv);
	
	int max_age = max(age);
	dvector iage(1,max_age);
	iage.fill_seqadd(1,1);
	dmatrix p_fl(1,nsex,1,max_age);
	for(int j=1;j<=nsex;j++)
	{
		for(int i=1;i<=max_age;i++)
		{
			dvariable l1  = exp(log_mu_l1(j));  // 50;       //
			dvariable l2  = exp(log_mu_l2(j));  // 150;      //
			dvariable rho = exp(log_mu_rho(j)); // exp(-0.2);//
			dvariable   b = exp(log_mu_b(j));   // 1.0;      //
			p_fl(j,i) = value(get_len(i,age_1,age_2,l1,l2,rho,b));
		}
	}
	REPORT(iage);
	REPORT(p_fl);
	
	//Transformation to standard von B parameters where
	// l_j = linf*(1-exp(-vbk*(age-to)))^p
	// p = 1/b
	// k = -log(rho)
	// linf = [(exp(vbk*t2)*l2^b-exp(vbk*t1)*l1^b))/(exp(vbk*t2)-exp(vbk*t1))]^(1/b)
	// to = t1 + t2 - 1/vbk*log[(exp(vbk*t2)*l2^b-exp(vbk*t1)*l1^b))/(exp(vbk*t2)-exp(vbk*t1))]
	
	dvector linf(1,nsex);
	dvector  vbk(1,nsex);
	dvector   t0(1,nsex);
	dvector    p(1,nsex);
	p    = 1./b;
	vbk  = -log(rho);
	dvector n1 = elem_prod(exp(vbk*age_2),pow(l2,b))-elem_prod(exp(vbk*age_1),pow(l1,b));
	dvector n2 = exp(vbk*age_2)-exp(vbk*age_1);
	linf = pow(elem_div(n1,n2),p);
	t0   = double(age_1)+double(age_2)-elem_prod(1.0/vbk,log(elem_div(n1,n2)));

	
	REPORT(linf);
	REPORT(vbk);
	REPORT(t0);
	REPORT(p);
	
	
FUNCTION mcmcReport
  {
	/*Write MCMC report*/
	static int nf=0;  /**< Counter for the number of function calls */
	
	if(nf==0)
	{
		cout<<"Writing MCMC report... Please wait"<<endl;
		ofstream ofs("vonBH.mcmc");
		ofs<<"sex \t area \t F0.1 \t M"<<endl;
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
			ofs<<j<<"\t"<<k<<"\t"<<f01(j,k)<<"\t"<<m<<endl;
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
	
	/*loop over sexes and areas and calculate F0.1*/
	double aa = 6.82e-6;
	double bb = 3.24;
	dvector m_m(1,2);
	m_m(1) = value(m);	//female m
	m_m(2) = 0.135/0.15*value(m);	//male m
	
	dvector la(1,nage);
	dvector wa(1,nage);
	dvector sa(1,nage);
	for(j=1;j<=nsex;j++)
	{
		for(k=1;k<=narea;k++)
		{
			/*Growth*/
			//double m_linf = value(linf(j,k));
			//double m_vbk  = value(vbk(j,k) );
			//double m_to   = value(to(j,k)  );
			//la = m_linf*(1.0-exp(-m_vbk*(age-m_to)));
			
	 	 	double m_l1  = value(l1(j,k));  // 50;       //
	 	 	double m_l2  = value(l2(j,k));  // 150;      //
	 	 	double m_rho = value(rho(j,k)); // exp(-0.2);//
	 	 	double   m_b = value(b(j,k));   // 1.0;      //
	 	 	for(i=1;i<=nage;i++)
	 	 	{
	 	 		la(i) =(get_len(i,age_1,age_2,m_l1,m_l2,m_rho,m_b));
	 	 	}
	 	 	wa = aa * pow(la,bb);
	 	 
	 	 	/*Selectivity*/
	 	 	sa = c_selexfun.linapprox(bin,CSelL(j),la);
	 	 
	 	 	/*F0.1 reference point*/
	 	 	yieldPerRecruit c_ypr(m_m(j),wa,sa);
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
  //  int i,j;
  //  int nage = max(age);
  //  dvector age(1,nage);
  //  age.fill_seqadd(1,1);
  //  dvector m(1,2);
  //  m(1) = 0.35;	//female m
  //  m(2) = 0.135;	//male m
  //  
  //  double a = 6.82e-6;
  //  double b = 3.24;
  //  dmatrix la(1,2,1,nage);
  //  dmatrix wa(1,2,1,nage);
  //  dmatrix lx(1,2,1,nage);
  //  dmatrix lz(1,2,1,nage);
  //  dmatrix sx(1,2,1,nage);
  //  Selex cselex;
  //  dvector bin(1,8);
  //  bin.fill("{50,60,70,80,90,100,110,120}");
  //  dmatrix CSelL(1,2,1,8);
  //  CSelL(1).fill("{1.63180e-09,3.25740e-09,0.0603022,0.300891,0.630344,0.913893,1,1}");
  //  CSelL(2).fill("{2.17159e-09,3.96878e-03,0.0567109,0.281436,0.585461,0.835614,1,1}");
  //  
  //  /*Growth*/
  //  la.initialize();
  //  wa.initialize();
  //  for(i=1;i<=2;i++)
  //  {
  //  	la(i) = value(linf(i,id)*(1.0-exp(-vbk(i,id)*(age-to(i,id)))));
  //  	wa(i) = a * pow(la(i),b);
  //  	sx(i) = cselex.linapprox(bin,CSelL(i),la(i));
  //  }
  //  
  //  //	yieldPerRecruit(double fe, double m, dvector wa, dvector sa);
  //  yieldPerRecruit c_ypr(m(1),wa(1),sx(1));
  //  //cout<<"YPR="<<c_ypr.getypr(0.01)<<endl;
  //  //cout<<"YPR="<<c_ypr.getypr(0.2)<<endl;
  //  //cout<<"YPR="<<c_ypr.getypr(0.3)<<endl;
  //  
  //  
  //  /*Survivorship*/
  //  lx.initialize();
  //  lz.initialize();
  //  lx = 1.0;
  //  lz = 1.0;
  //  for(i=1;i<=2;i++)
  //  {
  //  	for(j=2;j<=nage;j++)
  //  	{
  //  		lx(i,j) = lx(i,j-1)*exp(-m(i));
  //  		lz(i,j) = lz(i,j-1)*exp(-m(i)-fe*sx(i,j));
  //  	}
  //  }
  //  double ypr0 = c_ypr.getypr(0);
  //  //double ypr  = fe*sum(elem_prod(lz,elem_prod(wa,sx)));
  //  double ypr = c_ypr.getypr(fe);
  //  double f01 = c_ypr.getf01();
  //  cout<<fe<<"\t"<<ypr<<"\t"<<ypr0<<"\t"<<f01<<endl;
  //  return(ypr);
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
	

