// Class for Yield per recruit calculations and F0.1 and Fspr-based reference points.
#include <admodel.h>

#ifndef _YPR_H_
#define _YPR_H_
#define MAXITER 200
#define TOL 1.e-6

class yieldPerRecruit{
	
private:
	int    m_nage;
	int    m_sage;
	double m_fe;
	double m_f01;
	double m_m;
	double m_ypr;
	double m_dphiq;
	double m_d2phiq;
	double m_dphiq0;
	
	dvector m_lx;
	dvector m_lz;
	dvector m_dlz;
	dvector m_d2lz;
	dvector m_oa;
	dvector m_sa;
	dvector m_za;
	dvector m_qa;
	dvector m_wa;
	dvector m_va;
	
	
public:
	// default constructor
	yieldPerRecruit(double m, dvector wa, dvector va);
	~yieldPerRecruit(){} 	// descructor
	
	// getters
	double    getypr() { return m_ypr;     }
	double    getypr(double fe) { calcYieldPerRecruit(fe); return m_ypr; }
	double    getf01()          { calcF01(); return m_f01;}
	
	
	// member functions
	void calcSurvivorship(double fe, double m, dvector va);
	void calcYieldPerRecruit(double fe);
	void calcF01();
};

#endif



yieldPerRecruit::yieldPerRecruit(double m, dvector wa, dvector va)
{
	// Default constructor (one sex version)
	m_sage = wa.indexmin();
	m_nage = wa.indexmax();
	
	m_m  = m;
	m_wa = wa;
	m_va = va;
	
	calcYieldPerRecruit(0);
	m_dphiq0 = m_dphiq;		//initial slope at fe=0;
	//calcF01();
	
}

void yieldPerRecruit::calcF01()
{
	/*
		Calculate F0.1, which is defined as
		as the value of F where the tangential line
		on the yield per recruit plot is equal to 0.1 * dphiq0
		
		Use Newtons root finding method find
		fi = fi - (0.1*dphiq0-dphiq)/d2phiq
	*/
	
	double fi = 0.6*m_m;  // initial guess for F0.1
	int iter;
	double x1 = 0;
	double x2 = 10;
	
	for(iter=1;iter<=MAXITER;iter++)
	{
		calcYieldPerRecruit(fi);
		fi -= (0.1*m_dphiq0 - m_dphiq)/(m_d2phiq);
		// Check boundary conditions
		if((x1-fi)*(fi-x2)<0)
		{
			fi += 0.98*(0.1*m_dphiq0 - m_dphiq)/(m_d2phiq);
		}
		//cout<<"iter: "<<iter<<" fi="<<fi<<" dphiq="<<(0.1*m_dphiq0 - m_dphiq)<<endl;
		if(fabs(0.1*m_dphiq0 - m_dphiq) <= TOL) break;
	}
	
	m_f01 = fi;
	
}


void yieldPerRecruit::calcYieldPerRecruit(double fe)
{
	calcSurvivorship(fe,m_m,m_va);
	double phiq;
	double dphiq;
	double d2phiq;
	dvector qa(m_sage,m_nage);
	dvector t1(m_sage,m_nage);
	dvector t3(m_sage,m_nage);
	
	qa     = elem_div(elem_prod(elem_prod(m_wa,m_va),m_oa),m_za);
	phiq   = m_lz * qa;
	
	// 1st derivative for the yield per recruit
	t1     = elem_div(elem_prod(elem_prod(m_wa,m_va),m_sa),m_za);
	t3     = elem_div(elem_prod(qa,m_va),m_za);
	dphiq  = qa*m_dlz + m_lz*t1 + m_lz*t3;
	
	// 2nd derivative for per recruit yield (nasty)
	dvector t0 = elem_div(m_oa,m_za);
	dvector t2  = 2. * m_dlz;
	dvector V2  = elem_prod(m_va,m_va);
	dvector t5  = elem_div(elem_prod(m_wa,V2),m_za);
	dvector t7  = elem_div(m_va,m_za);
	dvector t9  = elem_prod(t5,m_sa);
	dvector t11 = elem_prod(t5,t0);
	dvector t13 = elem_prod(m_lz,t5);
	dvector t14 = elem_prod(m_va,m_sa);
	dvector t15 = elem_prod(t7,m_sa);
	dvector t17 = elem_prod(m_va,t0);
	dvector t18 = elem_div(t17,m_za);
	d2phiq      =  m_d2lz*qa + t2*t9 - t2*t11 - t13*t14 -2.*t13*t15 + 2.*t13*t18;


	
	// Member variables
	m_ypr    = fe*phiq;
	m_qa     = qa;
	m_dphiq  = dphiq;
	m_d2phiq = d2phiq;
}

void yieldPerRecruit::calcSurvivorship(double fe, double m, dvector va)
{
	int i;
	dvector   lx(m_sage,m_nage);
	dvector   lz(m_sage,m_nage);
	dvector  dlz(m_sage,m_nage);
	dvector d2lz(m_sage,m_nage);
	dvector   sa(m_sage,m_nage);
	dvector   za(m_sage,m_nage);
	dvector   oa(m_sage,m_nage);
	  lx(m_sage) = 1.0;
	  lz(m_sage) = 1.0;
	 dlz(m_sage) = 0.0;
	d2lz(m_sage) = 0.0;
	
	sa = exp(-m-fe*va);
	za = m+fe*va;
	oa = (1.-sa);
	for(i=m_sage+1; i<=m_nage; i++)
	{
		lx(i)   = lx(i-1)*exp(-m);
		lz(i)   = lz(i-1)*sa(i-1);
		dlz(i)  = sa(i-1)*(dlz(i-1)-lz(i-1)*va(i-1));
		d2lz(i) = sa(i-1) * ( d2lz(i-1)+ lz(i-1)*va(i-1)*va(i-1) );
		if(i==m_nage)
		{
			lx(i)  /= (1.-exp(-m));
			lz(i)  /= oa(i);
			dlz(i)  = dlz(i)/oa(i) - lz(i-1)*sa(i-1)*va(i)*sa(i)/square(oa(i));
			
			double V1  = va(i-1);
			double V2  = va(i);
			double oa2 = oa(i)*oa(i);
			
			d2lz(i) = d2lz(i)/oa(i) 
						+ 2*lz(i-1)*V1*sa(i-1)*V2*sa(i)/oa2
						+ 2*lz(i-1)*sa(i-1)*V2*V2*sa(i)*sa(i)/(oa(i)*oa2)
						+ lz(i-1)*sa(i-1)*V2*V2*sa(i)/oa2;
			
		}
	}
	m_lx   = lx;
	m_lz   = lz;
	m_dlz  = dlz;
	m_d2lz = d2lz;
	m_za   = za;
	m_sa   = sa;
	m_oa   = oa;
	
}