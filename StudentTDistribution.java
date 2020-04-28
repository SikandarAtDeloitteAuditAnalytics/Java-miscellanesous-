import javax.swing.JOptionPane;

public class StudentTDistribution {
	
	private double degreesOfFreedom = 2.0;

	public StudentTDistribution(double nu) {
		this.setDegreesOfFreedom(nu);
	}
	
	public void setDegreesOfFreedom(double nue) {
		if (nue < 0.0) {
			String txt = "t-distribution: Degrees of freedom erroneously set to value <= 0";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}
		this.degreesOfFreedom = nue;
	}

	public double getPDF(double x) {
		return calculateProbabilityDensityFunction(x);
	}

	public double getCDF(double x) {
		return calculateCumulativeDistributionFunction(x);
	}

	public double getInverseCDF(double p) {
		return calculateInverseCumulativeDistributionFunction(p);
	}

	
	
	/**
	 * Method calculateCumulativeDistributionFunction
	 * <p>
	 * calculates the cumulative distribution function (CDF) of Student's t distribution 
	 * 
	 * 
	 * @param x
	 *            = the point at which the pdf is to be evaluated
	 * 
	 * @return out
	 * 
	 * 
	 */
	private double calculateCumulativeDistributionFunction(double x) {
		double nu = this.degreesOfFreedom;
		double out = 0.0;
	
	            
		     if ((nu >= 2.0) || (nu < 1.0)) {
		       if (nu <= 100.0) {	 
		    	 
		    	 double arg = (x + Math.sqrt(x * x + nu)) / (2.0 * Math.sqrt(x * x + nu));
		    	 double prob = incompleteBetaFunction(nu/2.0, nu/2.0, arg);
                 out = prob;
                 
			     if (x == 0.0) {
				     out = 0.5;
			     }
			     if (x > 0.0) {
				     out = 1 - prob;
			     }
		       } else { 
		    	/*   out = cdft_aux(x); */
		    	  
		    	      double sd = Math.sqrt(nu/(nu-2.0));
		 			  double z = x / sd;
		 			  if( (nu < 10.0) && (z < -3000.0) ) { out = 0.0; return out; }
		 		      if( (nu >= 10.0) && (z < -150.0) ) { out = 0.0; return out; }
		 		      if( (nu < 10.0)  && (z > 3000.0) ) { out = 1.0; return out; }
		 		      if( (nu >= 10.0)  && (z > 150.0) ) { out = 1.0; return out; }
		 		      
		 		  
		 		  
		 		  if (nu > 100.0) {
		 			  double dcdfn = cumulativeDistributionFunctionStandardNormalScalar(x);
		 			  double d1 = x;
		 		      double d3 = Math.pow(x, 3.0);
		 		      double d5 = Math.pow(x, 5.0);
		 		      double d7 = Math.pow(x, 7.0);
		 		      double d9 = Math.pow(x, 9.0);
		 		      double d11 = Math.pow(x, 11.0);
		 		      double b11 = 0.25;
		 		      double b21 = 0.01041666666667; 
		 	          double b22 = 3.0; 
		 	          double b23 = -7.0;
		 	          double b24 = -5.0; 
		 	          double b25 = -3.0;
		 	 	      double b31 = 0.00260416666667;
		 	 	      double b32 = 1.0;
		 	 	      double b33 = -11.0;
		 	 	      double b34 = 14.0;
		 	 	      double b35 = 6.0;
		 	 	      double b36 = -3.0;
		 	 	      double b37 = -15.0;
		 	 	      double dconst = 0.3989422804; 
		 		      
		 		      double term1 = b11*(d3+d1)/nu;
		 		      double term2 = b21*(b22*d7+b23*d5+b24*d3+b25*d1)/(nu*nu);
		 		      double term3 = b31*(b32*d11+b33*d9+b34*d7+b35*d5+b36*d3+b37*d1)/(nu*nu*nu); 
		 		      double dcdf = term1+term2+term3;
		 		      dcdf = dcdfn - (dconst*(Math.exp(-x*x/2.0)))*dcdf;
		 		      out = dcdf;
		 		  }
		       }  
		     }
		     
		     if ((nu >= 1.0) && (nu < 2.0)) {
			     out = cdft_aux(x);
		     }
		
		
		return out;
	}

	
	/**
	 * Method calculateProbabilityDensityFunction
	 * 
	 * <p>
	 * calculates the probability density function (PDF) of Student's t distribution 
	 * 
	 * 
	 * @param x
	 *            = the point at which the pdf is to be evaluated
	 * 
	 * @return out
	 * 
	 * 
	 */
	
	private double calculateProbabilityDensityFunction(double x) {
		double nu = this.degreesOfFreedom;

		double out = gamma((nu + 1) / 2.0) / (Math.sqrt(Math.PI * nu) * gamma(nu / 2.0));
		out = out * Math.pow((1 + (x * x) / nu), -(nu + 1) / 2.0);
        
		return out;
	}

	
	/**
	 * Method cdft_aux
	 * 
	 * <p>
	 * calculation of the CDF of Student's t distribution via numerical integration
	 * using trapezoidal rule 
	 * 
	 * Source:https://en.wikipedia.org/wiki/Numerical_integration; see, in particular, the section
	 *        "Integrals over infinite intervals"
	 * 
	 * @param x
	 *            = the point at which the cdf is to be evaluated
	 * 
	 * @return out
	 * 
	 * 
	 */
	private double cdft_aux(double x) {
		double nu = this.degreesOfFreedom;
		double out = 0;
		double auxv1 = 0.00005;
		while (auxv1 < 1.0) {
			double b = Math.max(auxv1 + 0.00005, 1e-8);
			double a = Math.max(auxv1 - 0.00005, 1e-8);

			double ftb = 0.0;
			double fta = 0.0;

			double bb = x - ((1 - b) / b);
			ftb = gamma((nu + 1) / 2.0) / (Math.sqrt(Math.PI * nu) * gamma(nu / 2.0));
			ftb = ftb * Math.pow((1 + (bb * bb) / nu), -(nu + 1) / 2.0);
			ftb = ftb / (b * b);

			double aa = x - ((1 - a) / a);
			fta = gamma((nu + 1) / 2.0) / (Math.sqrt(Math.PI * nu) * gamma(nu / 2.0));
			fta = fta * Math.pow((1 + (aa * aa) / nu), -(nu + 1) / 2.0);
			fta = fta / (a * a);

			auxv1 = auxv1 + 0.0001;
			out = out + 0.5 * (b - a) * (ftb + fta);
		}
		return out;
	}


	/**
	 * Method calculateInverseCumulativeDistributionFunction
	 * 
	 * <p>
	 * calculation of the inverse CDF of Student's t distribution
	 * 
	 * 
	 * @param x
	 *            = the point at which the inverse of the cdf is to be evaluated
	 * 
	 * @return out
	 * 
	 * 
	 */
	
	private double calculateInverseCumulativeDistributionFunction(double p) {
		double nu = this.degreesOfFreedom;
		double out = 0.5;
		if (p <= 0.0) {
			String txt = "inverse CDF of t-distribution: p erroneously set to value <= 0";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}
		if (p >= 1.0) {
			String txt = "inverse CDF of t-distribution: p erroneously set to value >= 1";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}

		
		
		if ( ((nu >= 2.0) || (nu < 1.0)) && (nu <= 100.0) ) {
			double x = inverseOfIncompleteBetaFunction(2.0 * Math.min(p, 1.0 - p), 0.5 * nu, 0.5);
			x = Math.sqrt(nu * (1.0 - x) / x);
			out = x;
			if (p < 0.5) {
				out = -x;
			}
		}

		if ((nu >= 1.0) && (nu < 2.0)) {
			out = cdfti_aux(p);
		}
		
		if (nu > 100.0) {
			out = cdfti_aux(p);
		}
		
		return out;
	}

	/**
	 * Method cdfti_aux
	 * 
	 * <p>
	 * numerical, approximate calculation of the inverse CDF of Student's t distribution
	 * 
	 * 
	 * @param u
	 *            = the point at which the inverse of the cdf is to be evaluated
	 * 
	 * @return out
	 * 
	 * 
	 */
	
	private double cdfti_aux(double u) {
		double nu = this.degreesOfFreedom;
		
		if (u <= 0.0) {
			String txt = "Function cdfti_aux: u erroneously set to value <= 0";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}
		if (u >= 1.0) {
			String txt = "Function cdfti_aux: u erroneously set to value >= 1";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}
		
		
		int j = 0;
		double out = 0.0;
		int maxiter = 10000;
		double tol = 1e-8;
		double diff = 1e300;
		double trialx = 0.0;
		boolean done = false;
		while (done == false) {
			double fx = Math.abs(calculateCumulativeDistributionFunction(trialx) - u);
			diff = fx;
			if (diff < tol) {
				done = true;
				out = trialx;
			}

			if (j > maxiter && done == false) {
				String txt = "Calculation of inverse of Student t cdf at " + Double.toString(u) + " failed!";
				issueErrorMessageAndTerminateProgramExcecution(txt);
			}
			if (done == false) {
				double dfdx = Math.abs(calculateCumulativeDistributionFunction(trialx + 1e-8) - u);
				dfdx = dfdx - Math.abs(calculateCumulativeDistributionFunction(trialx - 1e-8) - u);
				dfdx = dfdx / (2.0 * 1e-8);
				trialx = trialx - fx / dfdx;
			}
			j = j + 1;
		}
		return out;
	}

	/**
	 * Method logGamma
	 * 
	 * <p>
	 * returns the log-Gamma function evaluated at x Calculation uses Lanczos
	 * approximation formula. Source:
	 * https://github.com/massiccio/java/blob/master/src/math/GammaFunction.java as
	 * of Sep 26, 2018
	 * 
	 * 
	 * @param x
	 *            = the point at which the Gamma function is to be evaluated
	 * 
	 * @return ln[Gamma(x)]
	 * 
	 * 
	 */

	private static double logGamma(double x) {

		double out = 0.0;

		if (x <= 0.0) {
			String txt = "Argument of logGamma function not > 0!";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}

		boolean xIsInteger = false;
		if (x == Math.floor(x)) {
			xIsInteger = true;
			if (x == 1.0) {
				out = 1.0;
			}
			if (x > 1.0) {
				out = 1.0;
				int j = 1;
				while (j < ((int) x)) {
					out = out * ((double) j);
					j = j + 1;
				}
			}
			out = Math.log(out);
		}

		if (xIsInteger == false) {
			// Log of Gamma from Lanczos with g=5, n=6/7
			// not in A & S
			double[] coef = { 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
					0.1208650973866179E-2, -0.5395239384953E-5 };

			double LogSqrtTwoPi = 0.91893853320467274178;

			double denom = x + 1;

			double y = x + 5.5;

			double series = 1.000000000190015;

			for (int i = 0; i < 6; ++i) {
				series += coef[i] / denom;
				denom += 1.0;
			}
			out = (LogSqrtTwoPi + (x + 0.5) * Math.log(y) - y + Math.log(series / x));

		}

		return out;
	}

	/**
	 * Method gamma
	 * 
	 * <p>
	 * returns the Gamma function evaluated at x
	 * 
	 * 
	 * @param x
	 *            = the point at which the Gamma function is to be evaluated
	 * 
	 * @return Gamma(x)
	 * 
	 */
	private static double gamma(double x) {
		return Math.exp(logGamma(x));
	}

	double hypergeometric2F1(double a, double b, double c, double x) {
		double TOLERANCE = 1.0e-10;
		double term = a * b * x / c;
		double value = 1.0 + term;
		double n = 1.0;

		while (Math.abs(term) > TOLERANCE) {
			a = a + 1.0;
			b = b + 1.0;
			c = c + 1.0;
			n = n + 1.0;
			term = term * (a * b * x / c / n);
			value = value + term;

		}

		return value;
	}

	/**
	 * Method incomplete Beta function
	 * 
	 * <p>
	 * returns the incomplete Beta function, evaluated at x and having parmeters a and b. 
	 * Source: Incomplete Beta function in C 
	 * https://codeplea.com/incomplete-beta-function-c
	 * of Apr 28, 2020
	 * 
	 * @param a
	 *            = first parameter
	 *            
	 * param b 
	 *            = second parameter           
	 * 
	 * @param x
	 *            = the point at which the incomplete beta function is to be evaluated
	 * 
	 * 
	 * 
	 */
	
	private static double incompleteBetaFunction(double a, double b, double x) {
		

		if (x < 0.0) {
			String txt = "Incorrect input (x < 0) to incomplete Beta function";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}
        
		if (x > 1.0) {
			String txt = "Incorrect input (x > 1) to incomplete Beta function";
			issueErrorMessageAndTerminateProgramExcecution(txt);
		}
		
		
		double  STOP = 1.0e-8;
		double  TINY = 1.0e-30;
		
	    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
	    if (x > (a+1.0)/(a+b+2.0)) {
	        return (1.0-incompleteBetaFunction(b,a,1.0-x)); /*Use the fact that beta is symmetrical.*/
	    }

	    /*Find the first part before the continued fraction.*/
	    double lbeta_ab = logGamma(a)+logGamma(b)-logGamma(a+b);
	    double front = Math.exp(Math.log(x)*a+Math.log(1.0-x)*b-lbeta_ab) / a;

	    /*Use Lentz's algorithm to evaluate the continued fraction.*/
	    double f = 1.0, c = 1.0, d = 0.0;

	    int i, m;
	    for (i = 0; i <= 10000; ++i) {
	        m = i/2;

	        double numerator;
	        if (i == 0) {
	            numerator = 1.0; /*First numerator is 1.0.*/
	        } else if (i % 2 == 0) {
	            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
	        } else {
	            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
	        }

	        /*Do an iteration of Lentz's algorithm.*/
	        d = 1.0 + numerator * d;
	        if (Math.abs(d) < TINY) d = TINY;
	        d = 1.0 / d;

	        c = 1.0 + numerator / c;
	        if (Math.abs(c) < TINY) c = TINY;

	        double cd = c*d;
	        f *= cd;

	        /*Check for stop.*/
	        if (Math.abs(1.0-cd) < STOP) {
	            return front * (f-1.0);
	        }
	    }

	    String txt = "Incomplete beta function: Maximum number of iterations (10 000) exceeded!";
	    issueErrorMessageAndTerminateProgramExcecution(txt);
	    return 1.0/0.0; /*Needed more loops, did not converge.*/

		
	}


	/**
	 * Method inverseOfIncomplete Beta function
	 * 
	 * <p>
	 * returns the inverse of the incomplete Beta function, evaluated at p and having parmeters a and b. 
	 * 
	 * @param a
	 *            = first parameter
	 *            
	 * param b 
	 *            = second parameter           
	 * 
	 * @param x
	 *            = the point at which the incomplete beta function is to be evaluated
	 * 
	 * 
	 * 
	 */
	
	private static double inverseOfIncompleteBetaFunction(double p, double a, double b) {
		double EPS = 0.000000001;
		double pp, t, u, err, x, al, h, w, afac;
		double a1 = a - 1.0;
		double b1 = b - 1.0;
		int j;
		x = -999.99;
		if (p <= 0.) {
			return 0.0;
		}
		if (p >= 1.) {
			return 1.0;
		}

		if ((p > 0.0) && (p < 1.0)) {
			if ((a >= 1.0) && (b >= 1.0)) {
				if (p < 0.5) {
					pp = p;
				} else {
					pp = 1.0 - p;
				}
				t = Math.sqrt(-2. * Math.log(pp));
				x = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
				if (p < 0.5) {
					x = -x;
				}
				al = ((x * x) - 3.0) / 6.0;
				h = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
				w = (x * Math.sqrt(al + h) / h)
						- (1. / (2. * b - 1) - 1. / (2. * a - 1.)) * (al + 5. / 6. - 2. / (3. * h));
				x = a / (a + b * Math.exp(2.0 * w));
			} else {
				double lna = Math.log(a / (a + b));
				double lnb = Math.log(b / (a + b));
				t = Math.exp(a * lna) / a;
				u = Math.exp(b * lnb) / b;
				w = t + u;
				if (p < t / w) {
					x = Math.pow(a * w * p, 1.0 / a);
				} else {
					x = 1. - Math.pow(b * w * (1.0 - p), 1.0 / b);
				}
			}
		}
		afac = -logGamma(a) - logGamma(b) + logGamma(a + b);

		j = 0;
		while (j < 10000) {
			if (x == 0.0 || x == 1.0) {
				return x;
			}
			err = incompleteBetaFunction(a, b, x) - p;
			t = Math.exp(a1 * Math.log(x) + b1 * Math.log(1.0 - x) + afac);
			u = err / t;
			t = u / (1.0 - 0.5 * Math.min(1.0, u * (a1 / x - b1 / (1.0 - x))));
			x = x - t;

			if (x <= 0.0)
				x = 0.5 * (x + t);
			if (x >= 1.0)
				x = 0.5 * (x + t + 1.0);
			if (Math.abs(t) < EPS * x) {
				if (j > 0) {
					j = 10000;
				}
			}
			j = j + 1;
		}
		return x;
	}

	/**
	 * Method cumulativeDistributionfunctionStandardNormalScalar
	 * 
	 * <p>
	 * returns the cumulative Standard Normal distribution valuated as z
	 * 
	 * 
	 * @param z
	 *            = the point at which the function is to be evaluated
	 * 
	 * @return cdfn = the resulting scalar
	 * 
	 */

	private static double cumulativeDistributionFunctionStandardNormalScalar(double z) {

		double a[][] = new double[5][1];
		double b[][] = new double[4][1];
		double c[][] = new double[9][1];
		double d[][] = new double[8][1];
		double p[][] = new double[6][1];
		double q[][] = new double[5][1];
		double xden;
		double y;
		double ysq;
		double xnum;
		double result = 0.00;
		double x;
		int i;

		double sqrpi = 0.564189583547756;
		double thresh = 0.46875;
		double sixten = 16;

		double xinf = 1.79E+308;
		double xneg = -26.628;
		double xsmall = 1.11E-16;
		double xbig = 26.543;
		double xhuge = 67100000;
		double xmax = 2.53E+307;
		double cdfn = -999.00;

		a[0][0] = 3.16112374387057;
		a[1][0] = 113.86415415105;
		a[2][0] = 377.485237685302;
		a[3][0] = 3209.37758913847;
		a[4][0] = 0.185777706184603;

		b[0][0] = 23.6012909523441;
		b[1][0] = 244.024637934444;
		b[2][0] = 1282.61652607737;
		b[3][0] = 2844.23683343917;

		c[0][0] = 0.56418849698867;
		c[1][0] = 8.88314979438838;
		c[2][0] = 66.1191906371416;
		c[3][0] = 298.6351381974;
		c[4][0] = 881.952221241769;
		c[5][0] = 1712.04761263407;
		c[6][0] = 2051.07837782607;
		c[7][0] = 1230.339354798;
		c[8][0] = 2.15311535474404E-08;

		d[0][0] = 15.7449261107098;
		d[1][0] = 117.693950891312;
		d[2][0] = 537.18110186201;
		d[3][0] = 1621.38957456669;
		d[4][0] = 3290.79923573346;
		d[5][0] = 4362.61909014325;
		d[6][0] = 3439.36767414372;
		d[7][0] = 1230.33935480375;

		p[0][0] = 0.305326634961232;
		p[1][0] = 0.360344899949804;
		p[2][0] = 0.125781726111229;
		p[3][0] = 1.60837851487423E-02;
		p[4][0] = 6.58749161529838E-04;
		p[5][0] = 1.63153871373021E-02;

		q[0][0] = 2.56852019228982;
		q[1][0] = 1.87295284992346;
		q[2][0] = 0.527905102951428;
		q[3][0] = 6.05183413124413E-02;
		q[4][0] = 2.33520497626869E-03;

		x = -z / 1.4142135623731;
		y = x;
		if (y < 0) {
			y = (-y);
		}

		if (y <= thresh) {
			ysq = 0;
			if (xsmall < y) {
				ysq = y * y;
			}
			xnum = a[4][0] * ysq;
			xden = ysq;
			i = 1;
			while (i <= 3) {
				xnum = (xnum + a[i - 1][0]) * ysq;
				xden = (xden + b[i - 1][0]) * ysq;
				i = i + 1;
			}
			result = 0.5 - 0.5 * x * (xnum + a[3][0]) / (xden + b[3][0]);

		}
		;
		if ((y > thresh) && (y <= 4)) {
			xnum = c[8][0] * y;
			xden = y;
			i = 1;
			while (i <= 7) {
				xnum = (xnum + c[i - 1][0]) * y;
				xden = (xden + d[i - 1][0]) * y;
				i = i + 1;
			}
			result = (xnum + c[7][0]) / (xden + d[7][0]);
			result = 0.5 * result * Math.exp(-(y * y));
			if (x < 0) {
				result = 1 - result;
			}
		}

		if ((y > 4)) {
			result = 0;
			if (y >= xbig) {
				cdfn = result;
				if (z > 0) {
					cdfn = 1;
				}
			}
			ysq = 1 / (y * y);
			xnum = p[5][0] * ysq;
			xden = ysq;
			i = 1;
			while (i <= 4) {
				xnum = (xnum + p[i - 1][0]) * ysq;
				xden = (xden + q[i - 1][0]) * ysq;
				i = i + 1;
			}
			result = ysq * (xnum + p[4][0]) / (xden + q[4][0]);
			result = (sqrpi - result) / y;
			result = 0.5 * result * Math.exp(-(y * y));
			if (x < 0) {
				result = 1 - result;
			}

		}
		if (cdfn == -999.00) {
			cdfn = result;
		}
		if (cdfn == -999.00) {
			System.out.println(z + " " + cdfn);
		}
		return cdfn;
	}


	private static void issueErrorMessageAndTerminateProgramExcecution(String txt) {

		JOptionPane.showOptionDialog(null, txt, "ERROR", JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,
				null, new Object[] {}, null);
		System.exit(0);
	}

}
