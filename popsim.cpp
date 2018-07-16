//[[Rcpp::depends(RcppArmadillo,RcppGSL)]]

#include <RcppArmadillo.h>
#include <omp.h>
#include<RcppGSL.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>



// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;


//Generated random deviate from truncated normal using inverse cdf method
//Adapted from John Burkardt code to armadillo syntax
//[[Rcpp::export]]
double rtnest( double mu, double s, double a, double b){
 
  double alpha = ( a - mu ) / s;
  double beta = ( b - mu ) / s;

  double alpha_cdf = gsl_cdf_ugaussian_P( alpha );
  double beta_cdf = gsl_cdf_ugaussian_P( beta );

  double u = as_scalar(randu(1));
  double xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf );
  double xi = gsl_cdf_ugaussian_Pinv( xi_cdf );
  
  if(is_finite(xi)){
	xi = xi;
  }
  else{
	xi =0;
  }
  
  double x = mu + s * xi;
  if(x<a){
	x = a;
  }
   if(x>b){
	x = b;
  }
  
  return x;  
}

//[[Rcpp::export]]
int rint(const vec prob, const ivec values){
	//Function to determine a random integer based on probability vector
	//Arguments:
		//prob - probability of each vector
		//values - vector of values to select from (note: same order as prob) 
	int out = 0;
	double test = as_scalar(randu(1));
	vec probuse = prob/sum(prob); //make sure prob. vector sums to 1
	vec cumprob = cumsum(probuse);
	uvec intv = find(cumprob <= test); 
	if(intv.n_elem>=1){
		out = max(intv)+1;
	}else{
		out =0;
	}	
	//Rcpp::Rcout<<cumprob<<":"<<intv<<":"<<test<<std::endl;
	return(out);
}

//[[Rcpp::export]]
List nxtyr(const vec Nmale, const vec Nfemale, const mat hazmale, const mat hazfemale,
			const vec rprob, const vec fawnprob, const double sexpr, const double report,
			const double saprob, const int maxfawn, const int years, const int cores){
	//Function to generate age-specific population and harvest sizes
	//Arguments:
		//Nmale - vector of current number of male individuals in each class of interest
		//Nfemale - vector of current number of female individuals in each class of interest
		//haz - matrix of hunting and non-hunting hazards for each class for each sex
		//rprob - age-specific birth probability
		//fawnprob - probability of the number of fawns for each doe [0:maxfawn]
		//sexpr - prob fawn is male
		//report - prob of reporting harvest
		//saprob - prob harvested animal sexed and aged
		//maxfawn - maximum number of fawns per doe
		//years - number of years to run process
		//cores - number of threads for parallel processing via openMP

	omp_set_num_threads(cores);  //determine number of cores

	mat Npremale(years,Nmale.n_elem,fill::zeros); //pre-hunt male population
	mat Nprefemale(years,Nfemale.n_elem,fill::zeros); //pre-hunt female population
	mat Npostmale(years,Nmale.n_elem,fill::zeros); //post-hunt male population at end of biological year (4/14)
	mat Npostfemale(years,Nfemale.n_elem,fill::zeros); //post-hunt female population at end of biological year (4/14)
	
	double mfawn=0;  //male fawns
	double ffawn=0;  //female fawns
	
	mat harvestm (years,Nmale.n_elem,fill::zeros); //harvested male population
	mat harvestf (years,Nfemale.n_elem,fill::zeros); //harvested female population
	mat oharvest (years,2,fill::zeros); //reported harvest - first column=antlered/ second column = antlerless
	
	mat saharvf(years,Nfemale.n_elem,fill::zeros); //reported and sex/age -females
	mat saharvm(years,Nmale.n_elem,fill::zeros); //reported and sex/age -males

for(int kk=0;kk<years;++kk)	{
	vec Nmaleuse(Npostmale.n_cols,fill::zeros);
	vec Nfemaleuse(Npostfemale.n_cols,fill::zeros);
	
	if(kk==0){
		Nmaleuse = Nmale;
		Nfemaleuse = Nfemale;
		Npremale(0,span::all) =  trans(Nmale); //set so initial values are prehunt population
		Nprefemale(0,span::all) = trans(Nfemale);
	}else{
		Nmaleuse(0) = 0; //no fawns to start year
		Nfemaleuse(0) = 0; //no fawns to start year
		Nmaleuse(span(1,Npostmale.n_cols-1)) = trans(Npostmale(kk-1,span(0,Npostmale.n_cols-2)));
		Nfemaleuse(span(1,Npostfemale.n_cols-1)) = trans(Npostfemale(kk-1,span(0,Npostfemale.n_cols-2)));
		Nmaleuse(Npostmale.n_cols-1) = Nmaleuse(Nmaleuse.n_elem-1)+Npostmale(kk-1,Npostmale.n_cols-1); //accumulator age class
		Nfemaleuse(Npostfemale.n_cols-1) = Nfemaleuse(Nfemaleuse.n_elem-1)+Npostfemale(kk-1,Npostfemale.n_cols-1); //accumulator age class
	}
		
	
	//Adults
	for(uword i=1; i<Nmaleuse.n_elem; ++i){

			mfawn=0; //reset to zero
			ffawn=0;
			vec test1 = (-1/(hazmale(0,i)))*log(randu(Nmaleuse(i))); //inverse-cdf method - exp random variable for non-hunting season
			vec test2 = (-1/(hazmale(1,i)))*log(randu(Nmaleuse(i))); //inverse-cdf method - exp random variable for hunting season
			vec test3 = (-1/(hazmale(0,i)))*log(randu(Nmaleuse(i))); //inverse-cdf method - exp random variable for non-hunting season
			
			vec test4 = (-1/(hazfemale(0,i)))*log(randu(Nfemaleuse(i))); //inverse-cdf method - exp random variable for non-hunting season
			vec test5 = (-1/(hazfemale(1,i)))*log(randu(Nfemaleuse(i))); //inverse-cdf method - exp random variable for hunting season
			vec test6 = (-1/(hazfemale(0,i)))*log(randu(Nfemaleuse(i))); //inverse-cdf method - exp random variable for non-hunting season
			
			mat bdate(Nfemaleuse(i),3,fill::zeros);  //birth date of fawns

			//Males
			#pragma omp for
			for(uword j = 0; j<test1.n_elem; ++j){
				if(kk>0){ //year 1 use initial values for prehunt populations
					if(test1(j) <= 169.0){  //Survive from 4/15 to 10/1 //Need 169.0 so returns double rather than integer (i.e., 0)
						continue; //skip to next individual because of death
					}else{
						Npremale(kk,i) = Npremale(kk,i)+1;
					}
				}
					
				if(test2(j) <= 92.0){ //Harvested during hunting season

						harvestm(kk,i) = harvestm(kk,i)+1;
						if(as_scalar(randu(1)) <= report){
							oharvest(kk,0) = oharvest(kk,0)+1;
							if(as_scalar(randu(1)) <= saprob){
								saharvm(kk,i) = saharvm(kk,i)+1;
							}
						}
						continue;
 				}
				
				if(test3(j) > 104.0){ //Survive from 1/1 to 4/15
					Npostmale(kk,i) = Npostmale(kk,i)+1;
				}
			}

			//Females
			#pragma omp for
			for(uword j=0; j<test4.n_elem; ++j){
				int numfawn =  rint(fawnprob,randi(1,distr_param(0,maxfawn))); //number of fawns for each doe
				if(as_scalar(randu(1)) <= rprob(i)){
					bdate(j,0) = rtnest(151,10,120,182); //fawn birth date (peak-5/31 and range from 4/30-7/1
				}else{
					bdate(j,0) = 10000000;
				}
				if((test4(j)+105) > bdate(j,0) && kk > 0){
					for(int f=0; f<numfawn; ++f){
						if(as_scalar(randu(1)) < sexpr){
							mfawn = mfawn+1;
							bdate(j,1) = bdate(j,1)+1; //indicator of male fawn
						}else{
							ffawn = ffawn+1;
							bdate(j,2) = bdate(j,2)+1; //indicator of male fawn
						}
					}
				}
				if(kk > 0){
					if(test4(j) <= 169.0){  
						continue;
					}else{ //Survive from 4/15 to 10/1
						Nprefemale(kk,i) = Nprefemale(kk,i)+1;
					}
				}
				if(test5(j) <= 92.0){ //Harvested during hunting season
						harvestf(kk,i) = harvestf(kk,i)+1;
						if(as_scalar(randu(1)) <= report){
							oharvest(kk,1) = oharvest(kk,1)+1;
							if(as_scalar(randu(1)) <= saprob){
								saharvf(kk,i) = saharvf(kk,i)+1;
							}
						}
						continue;
  					}
				if(test6(j) > 104.0){ //Survive from 1/1 to 4/15
						Npostfemale(kk,i) = Npostfemale(kk,i)+1;
					}
			}
			
			if(kk == 0 && i == 1){ //set to initial value for first year
				mfawn = Nmale(0);
				ffawn = Nfemale(0);
			}
			if(kk == 0 && i != 1){ //no fawns for other age classes, all captured in first year and are the initial values provided
				mfawn = 0;
				ffawn = 0;	
			}
		
			mat bdateextend(1,bdate.n_cols,fill::zeros);
			int lll=0;
			
			//Add birth dates for each fawn as rows
			if(kk >0){
				#pragma omp for	
				for(uword l=0; l<bdate.n_rows; ++l){
					if(bdate(l,1) >= 1){		
						for(int ll=0; ll<as_scalar(bdate(l,1)); ++ll){						
							lll++;
							rowvec tempmind(3);
							tempmind(0) = bdate(l,0);
							tempmind(1) = 1; //add male indicator only
							tempmind(2) =0;
							bdateextend.insert_rows(lll,tempmind);
						}
					}
					
					if(bdate(l,2) >= 1){				
						for(int ll=0; ll<as_scalar(bdate(l,2)); ++ll){						
							lll++;
							rowvec tempind(3);
							tempind(0) = bdate(l,0);
							tempind(1) = 0; //add female indicator only
							tempind(2) =1;
							bdateextend.insert_rows(lll,tempind);
						}
					}
					if(accu(bdate(l,span(1,2)))==0){
						lll++;
						bdateextend.insert_rows(lll, bdate(l, span::all));
					}	
				}
				bdateextend.shed_row(0); //remove initialization row
			}		
			//Male Fawns
			vec bdateusem = bdateextend.col(0);
			bdateusem = bdateusem(find(bdateextend(span::all,0) < 10000000 && bdateextend(span::all,1) != 0)); //subset birth dates for living fawns

			vec test7 = (-1/(hazmale(0,0)))*log(randu(mfawn)); //inverse-cdf method - exp random variable for non-hunting season
			vec test8 = (-1/(hazmale(1,0)))*log(randu(mfawn)); //inverse-cdf method - exp random variable for hunting season
			vec test9 = (-1/(hazmale(0,0)))*log(randu(mfawn)); //inverse-cdf method - exp random variable for non-hunting season			
			
			#pragma omp for
			for(uword j = 0; j<test7.n_elem; ++j){
				if(kk > 0){
					if(test7(j) <= (274-bdateusem(j))){  //Survive from birth to 10/1
						continue;
					}else{
						Npremale(kk,0) = Npremale(kk,0)+1;
					}
				}
				
				if(test8(j) <= 92.0){ //Harvested during hunting season
						harvestm(kk,0) = harvestm(kk,0)+1;
						if(as_scalar(randu(1)) <= report){
							
							
							
//Remove if only doing females because not including male fawns in antlerless count	//	oharvest(kk,1)=oharvest(kk,1)+1;  
							if(as_scalar(randu(1)) <= saprob){
								saharvm(kk,0) = saharvm(kk,0)+1;
							}
						}
						continue;
  				}
				if(test9(j) > 104.0){ //Survive from 1/1 to 4/15
								Npostmale(kk,0) = Npostmale(kk,0)+1;
							}
			}
		
			//Female Fawns
			vec bdateusef = bdateextend.col(0);
			bdateusef = bdateusef(find(bdateextend(span::all,0)<10000000 && bdateextend(span::all,2)!=0)); //subset birth dates for living fawns
			
			vec test10 = (-1/(hazfemale(0,0)))*log(1-randu(ffawn)); //inverse-cdf method - exp random variable for non-hunting season
			vec test11 = (-1/(hazfemale(1,0)))*log(1-randu(ffawn)); //inverse-cdf method - exp random variable for hunting season
			vec test12 = (-1/(hazfemale(0,0)))*log(1-randu(ffawn)); //inverse-cdf method - exp random variable for non-hunting season
			
			#pragma omp for	
			for(uword j=0; j<test10.n_elem; ++j){
				if(kk > 0){	
					if(test10(j) <= (274-bdateusef(j))){  //Survive from birth to 10/1
						continue;
					}else{
						Nprefemale(kk,0) = Nprefemale(kk,0)+1;
					}
				}
				if(test11(j) <= 92.0){ //Harvested during hunting season
					harvestf(kk,0) = harvestf(kk,0)+1;
					if(as_scalar(randu(1)) <= report){
						oharvest(kk,1) = oharvest(kk,1)+1;
						if(as_scalar(randu(1)) <= saprob){
							saharvf(kk,0) = saharvf(kk,0)+1;
						}
					}
						continue;
  					}
				if(test12(j) > 104.0){ //Survive from 1/1 to 4/15
						Npostfemale(kk,0) = Npostfemale(kk,0)+1;
					}
			}
	}
}
	return List::create(Named("Nprefemale") = Nprefemale,
                        Named("Npremale") = Npremale,
						Named("Npostmale") = Npostmale,
                        Named("Npostfemale") = Npostfemale,
						Named("harvestm") = harvestm,
						Named("harvestf") = harvestf,
						Named("oharvest") = oharvest,
						Named("saharvm") = saharvm,
						Named("saharvf") = saharvf);
		
}
