//Emily Morgan Liljestrand
//Releases and returns  divided into 4 zones, 12 months, 3.5 years
//Zone specific tagging loss/mortality based on "additional observed predation"
//6 movement matrices describe the movement on a monthly basis (5 matrices) and half-year basis (1 matrix)
//Priors placed on migration matrix values according to a dirchlet distribution assuming passive diffusion
//Priors on the catchability values according to normal distribution of the deviations of the mean
//A uniform prior distribution placed on the log scale monthly natural mortality rate
//The overdispersion parameter, k of the negative binomial distribution around recoveries, is fixed
//Created: May 3, 2017
//Updated: June 19, 2017

TOP_OF_MAIN_SECTION
  arrmblsize = 2000000; //sets the amount of memory used by the model

DATA_SECTION
  init_int numsts //number of monthly time steps
  init_int numsea //number of months
  init_int numzone //number of regions
  init_matrix releases(1,numsts,1,numzone) //release matrix
  init_ivector release_sea(1,numsts) //a vector with the month identity (1-January, 2- February, 3-March, etc.) for each column of "releases"
  init_ivector release_yr(1,numsts) //a vector with the year identity for each column of "releases" 1-66, 2-67, 3-68, 4-69
  init_vector tagsurv(1,numzone) //survival after tagging mortality/shedding for each region
  init_matrix effort(1,numsts,1,numzone)//effort (week*ton*vessels) for each time step and region
  init_matrix mageff(1,numsts,1,numzone)//vector of magnet efficiency for each month and region
  init_4darray recapture(1,numsts,1,numzone,1,numsts,1,numzone)//4d array of observed recaptures
  init_int effsamplesize //effective sample size used for the Dirchlet distribution
  init_3darray migrationprior(1,2,1,numzone,1,numzone) //prior values for the monthly and half-year movement matrices
  init_int numtheta //the number of theta parameters in this model. Should be (numzone*numsea)-(number of months and regions with no effort)-(numzone)
  init_int nummageff //the number of mageff parameters in this model
  init_vector mageff_bound(1,2) //4 regions median eff
  init_int testnum //test if the end of file reached

  int z //looping variable for region
  int s //looping variable for months releases
  int s2//looping variable for months recaptures
  int z2 //another looping variable for region
  int i//looping variable for all purposes

 LOCAL_CALCS
  if(testnum!=1234) //If the final value is not 1234, like it should be, the code stops and shows all the read-in values
  {
    cout<<"Num monthly time steps:"<<endl<<numsts<<endl;
    cout<<"Num months:"<<endl<<numsea<<endl;
    cout<<"Num zones:"<<endl<<numzone<<endl;
    cout<<"Release matrix:"<<endl<<releases<<endl;
    cout<<"Release season:"<<endl<<release_sea<<endl;
    cout<<"Release year:"<<endl<<release_yr<<endl;
    cout<<"Tagging survival:"<<endl<<tagsurv<<endl;
    cout<<"Effort Matrix:"<<endl<<effort<<endl;
    cout<<"Magnet Efficiency:"<<endl<<mageff<<endl;
    cout<<"Recapture Matrix:"<<endl<<recapture<<endl;
    cout<<"Effective Sample Size:"<<endl<<effsamplesize<<endl;
    cout<<"Migration Matrix Priors:"<<endl<<migrationprior<<endl;
    cout<<"Number of Theta parameters:"<<endl<<numtheta<<endl;
    cout<<"Test number:"<<endl<<testnum<<endl;
    exit(1);
  }
 END_CALCS

PARAMETER_SECTION
  init_bounded_number log_M(-5,5,1) //log natural mortality
  //number log_M
  init_bounded_vector log_q(1,numzone,-10,0,1) //log regional catchability
  init_vector theta(1,numtheta,1) //theta term for describing the month-by-region catchability
  init_bounded_number log_k(-5,3,-1) //log overdispersion paramter, (-1) indicates that it's not estimated
  init_3darray logmigparam(1,numsea-6,1,numzone-1,1,numzone,1) //log scale parameters for movement. The number of columns numzone-1 because one movement value is calculated by subtracting 1 from the others
  init_bounded_vector log_mageffv(1,nummageff,-3.5,-0.05,1) // log mageff

  //The non-log ajusted parameter values:
  number M //Monthly natural mortality
  number Myr //yearly natural mortality
  vector q(1,numzone) //regional catchability
  vector mageffv(1,nummageff) // mageff vector
  number k //k parameter in negative binomial

  //Total and month-by-region catchability
  matrix p(1,numsea,1,numzone) //monthly component of total catchability p=0 if there is no effort, p=1 in July, p=e^theta for all other values 
  matrix qtot(1,numsea,1,numzone)//matrix of total catchability (both monthly and zone components)

  matrix F(1,numsts,1,numzone) //fishing mortality
  matrix Z(1,numsts,1,numzone) //total mortality (F+M)
  matrix pobs(1,numsts,1,numzone) //probability of being observed in that region/time

  matrix releases_tm(1,numsts,1,numzone) //releases after tagging mortality applied
  matrix catches(1,numsts,1,numzone) //number estimated recapture each month for each zone after mageff
  matrix orecapture(1,numsts,1,numzone) //number obs recapture each month for each zone after mageff
  matrix N_pop(1,numsts,1,numzone) //number estimated tag_population each month for each zone
  matrix est_mageff(1,numsts,1,numzone) // estimated mean magnet efficiency

  3darray migmat(1,numsea-6,1,numzone,1,numzone)//the movement matrices, one for each month
  4darray Ntagmat(1,numsts,1,numzone,1,numsts,1,numzone)//matrix of matrices, for the number of tagged indiviuals still alive from each cohort at subsequent times

  number summ //a sum for adding all values in the column of migration matrix, for normalizing to <1 and sum to 1
  number summ0 //a sum 
  number count //a count variable for properly indicating which migration matrix to run

  4darray expected(1,numsts,1,numzone,1,numsts,1,numzone) //the expected number of tag returns based on movement and mortality parameters
  4darray likeli(1,numsts,1,numzone,1,numsts,1,numzone) //4d array of all likelihood values associated with each observation

  //These are necessary to get the mean catchability value when you aren't sure how many values there actually are
  number sum_qtot //summing all the total catchability (qtot values)
  number count_qtot //counting the number of total catchability (qtot values)
  number mean_qtot //calculating the mean based on the sum and the count

  number count_sdrepq //counter to effectively loop through the sd report values for q

  //These sdreport terms are necessary to understand the likelihood profile, including 95% C.I., of values that aren't actual parameters
  sdreport_matrix migmatsolutions1(1,numzone,1,numzone)//
  sdreport_matrix migmatsolutions2(1,numzone,1,numzone)
  sdreport_matrix migmatsolutions3(1,numzone,1,numzone)
  sdreport_matrix migmatsolutions4(1,numzone,1,numzone)
  sdreport_matrix migmatsolutions5(1,numzone,1,numzone)
  sdreport_matrix migmatsolutions6(1,numzone,1,numzone)
  sdreport_vector qtotsolutions(1,numtheta+numzone)

  //Likelihood components
  number negLLrecapture //Likelihood related to the number of recaptured tags
  number negLLmig //Penalty of the prior distribution on movement values
  number negLLq //Penalty of the prior distribution on deviation from the mean of total catchability values
  objective_function_value negLL //Objective function, the sum of all likelihood components

 LOCAL_CALCS
  //This section sets initial values for the parameters
  log_M=log(0.001);
  for(s=1;s<=numzone;s++)
  {
    log_q(s)=log(0.001);
  }
  log_k=log(2.5);
  for(s=1;s<=numsea-6;s++)
  {
    logmigparam(s)=1.0;
  }
  for(s=1;s<=numtheta;s++)
  {
    theta(s)=1;
  }

  for(s=1;s<=nummageff;s++)
  {
    log_mageffv(s)=log(0.4);
  }

  for(s=1;s<=numsea;s++)
  {
    for(z=1;z<=numzone;z++)
    {
    est_mageff(s,z)<=mageff(s,z);
    }
  }
 END_CALCS
  

PROCEDURE_SECTION
  exp_param(); //Exponentiates the log paramters
  calc_mageff();
  app_tagsurv(); //Apply the tagging survival rate to the observed returns 
  calc_qtot(); //Calculates total catchability from month and zone specific values
  calc_mort(); //Calculates natural, fishing, and total mortality
  calc_pobs(); //Recalculates probability of observation according to zone specific magnet efficiency
  calc_migmat(); //Calculates the migration matrices
  calc_ntagmat(); //Calculates the number tagged in each month and zone
  calc_negLL(); //Calculates the negLL based on observed and expected tag recoveries, migmat priors, and theta priors
  //calc_mageff();

  set_solutions(); //Calculates "migrationmatrixsolutions" which can be easily printed out to rep file
  
  /*
  cout<<negLLrecapture<<endl;
  cout<<negLLmig<<endl;
  cout<<negLLq<<endl;
  cout<<likeli<<endl;  
  exit(1);
  */

  if(mceval_phase())
  {
    write_waic_row();   // NEW: dump pointwise log-lik for this draw
    mcmcoutput();
  }
  
FUNCTION exp_param
  M=exp(log_M); //changes from log scale to actual M value
  for(s=1;s<=numzone;s++)
  {
    q(s)=exp(log_q(s)); //changes from log scale to actual q value
  }
  k=exp(log_k); //changes from log scale to actual k value
  
  //cout<<"M:"<<M<<endl;
  //cout<<"q:"<<q<<endl;
  //cout<<"k:"<<k<<endl;
  //cout<<"Numsea: "<<numsea<<endl;
  //exit(1);

FUNCTION app_tagsurv
  for(s=1;s<=numsts;s++)
  {
    releases_tm(s)=elem_prod(releases(s),tagsurv); //Applying tagging survival
  }
  
  //cout<<"Releases after tagging survival:"<<endl<<releases_tm<<endl;
  //exit(1);

FUNCTION calc_qtot
  //convert the theta values (parameters) into the p values (monthly catchability) and then multiply with zone catchability to get 4 x 12 matrix
  p(1,1)=0;
  p(2,1)=0;
  p(3,1)=0;
  p(4,1)=0;
  p(5,1)=exp(theta(1));
  p(6,1)=exp(theta(2));
  p(7,1)=1;
  p(8,1)=exp(theta(3));
  p(9,1)=exp(theta(4));
  p(10,1)=exp(theta(5));
  p(11,1)=exp(theta(6));
  p(12,1)=0;
  
  p(1,2)=0;
  p(2,2)=0;
  p(3,2)=0;
  p(4,2)=0;
  p(5,2)=exp(theta(7));
  p(6,2)=exp(theta(8));
  p(7,2)=1;
  p(8,2)=exp(theta(9));
  p(9,2)=exp(theta(10));
  p(10,2)=exp(theta(11));
  p(11,2)=exp(theta(12));
  p(12,2)=exp(theta(13));
  
  p(1,3)=exp(theta(14));
  p(2,3)=0;
  p(3,3)=0;
  p(4,3)=exp(theta(15));
  p(5,3)=exp(theta(16));
  p(6,3)=exp(theta(17));
  p(7,3)=1;
  p(8,3)=exp(theta(18));
  p(9,3)=exp(theta(19));
  p(10,3)=exp(theta(20));
  p(11,3)=exp(theta(21));
  p(12,3)=exp(theta(22));
  
  p(1,4)=0;
  p(2,4)=0;
  p(3,4)=0;
  p(4,4)=exp(theta(23));
  p(5,4)=exp(theta(24));
  p(6,4)=exp(theta(25));
  p(7,4)=1;
  p(8,4)=exp(theta(26));
  p(9,4)=exp(theta(27));
  p(10,4)=exp(theta(28));
  p(11,4)=exp(theta(29));
  p(12,4)=0;

  for(s=1;s<=numsea;s++)
  {
    qtot(s)=elem_prod(p(s),q); //Calculate total catchability from regional (q) effect, and month-by-region (p(s)) effect
  }
  //cout<<"Catchability Matrix: "<<endl<<qtot<<endl;
  //exit(1);

FUNCTION calc_mort
  for(s=1;s<=numsts;s++)
  {
    F(s)=elem_prod(effort(s),qtot(release_sea(s))); //calculate fishing mortality from effort and region/month specific catchability
  }
  
  //cout<<"F: "<<F<<endl;
  //exit(1);
  
  Z=F+M; //total mortality, the sum of fishing and natural mortality
  Myr=M*12; //yearly natural mortality, the monthly natural mortality times 12
  
  //cout<<"F:"<<endl<<F<<endl;
  //cout<<"Z:"<<endl<<Z<<endl;
  //cout<<"Effort:"<<endl<<effort<<endl;
  //exit(1);
  
  pobs=elem_prod(elem_div(F,Z),(1.0-exp(-Z))); //baranov catch equation, calculates the probability of CAPTURING a tagged individual based on total and fishing mortality

  //cout<<"Observation Before Mag Eff:"<<endl<<pobs<<endl;
  //exit(1);

FUNCTION calc_pobs
  for(s=1;s<=numsts;s++)
  {
    pobs(s)=elem_prod(pobs(s),est_mageff(s)); //adjusts the pobs even further by appling the region-specific magnet efficiency. This is the probability of CAPTURING and DETECTING a tagged individual
  }
  
  //cout<<"Total Probability of Observation:"<<endl<<pobs<<endl;
  //cout<<"Log Migration Parameters:"<<endl<<logmigparam<<endl;
  //exit(1);
  

FUNCTION calc_migmat
  //loop to calculate the regular movement matrix from the log-adjusted value. The non-estimated value that is calculated from the other values changes for ease of estimation/parameterizaion
  for(s=1;s<=numsea-6;s++)
  {
    for(z=1;z<=numzone;z++)
    {
      summ=0.0; //sets an initial sum for all values in a column, which ultimate must sum to 1
      if(z<=2) //if the column (region at the beginning of the month) is 1 or 2, set the movement values such that movement to region 2 is calculated, not estimated
        {
        for(z2=1;z2<=numzone;z2++) //for this matrix and column, we are going row by row
        {
          if(z2==1)
          {
            summ+=exp(logmigparam(s,z2,z)); //take each element in the parameter array, exponentiates it, and adds it to the growing sum
          }
          if(z2==2)
          {
            summ+=exp(1.0);
          }
          if(z2>2)
          {
            summ+=exp(logmigparam(s,z2-1,z));
          }
        //cout<<exp(logmigparam(s,z2,z))<<endl;
        }
        for(z2=1;z2<=numzone;z2++) //now that we have the sum, setting each value in the migration matrix by dividing each element by this sum
        {
          if(z2==1)
          {
            migmat(s,z2,z)=exp(logmigparam(s,z2,z))/summ; //for each element in the parameter array, take the term and divide by the total to get the proportion
          }
          if(z2==2)
          {
            migmat(s,z2,z)=exp(1.0)/summ;
          }
          if(z2>2)
          {
            migmat(s,z2,z)=exp(logmigparam(s,z2-1,z))/summ;
          }
        //cout<<migmat(s,z2,z)<<endl;
        }
      }
      
      if(z>2)//if the column (region at the beginning of the month) is 3 or 4, set the movement values such that movement to region 3 is calculated, not estimated
      {
        for(z2=1;z2<=numzone;z2++) //for this matrix and column, we are going row by row
        {
          if(z2<=2)
          {
            summ+=exp(logmigparam(s,z2,z)); //take each element in the parameter array, exponentiates it, and adds it to the growing sum
          }
          if(z2==3)
          {
            summ+=exp(1.0);
          }
          if(z2==4)
          {
            summ+=exp(logmigparam(s,z2-1,z));
          }
          //cout<<exp(logmigparam(s,z2,z))<<endl;
        }
        for(z2=1;z2<=numzone;z2++) //now that we have the sum, setting each value in the migration matrix by dividing each element by this sum
        {
          if(z2<=2)
          {
            migmat(s,z2,z)=exp(logmigparam(s,z2,z))/summ; //for each element in the parameter array, take the term and divide by the total to get the proportion
          }
          if(z2==3)
          {
            migmat(s,z2,z)=exp(1.0)/summ;
          }
          if(z2==4)
          {
            migmat(s,z2,z)=exp(logmigparam(s,z2-1,z))/summ;
          }
        //cout<<migmat(s,z2,z)<<endl;
        }
      }
      //migmat(s,numzone,z)=exp(1.0)/summ; //set the final value in the column to exp(1)/sum
    }
  }
  //cout<<"Migration Matrices:"<<endl<<migmat<<endl;
  //exit(1);

FUNCTION calc_ntagmat
  //calculates the number of survivng individuals in each cohort at each time after release
  for(s=1;s<=numsts;s++)
  {
    for(z=1;z<=numzone;z++)
    {
      Ntagmat(s,z)=0.0; //setting the initial number of tagged invididuals in subsequent years to 0
      Ntagmat(s,z,s,z)=releases_tm(s,z); //set the known number of releases at the initial region and month
      for(s2=s+1;s2<=numsts;s2++) //calculate the number of expected tagged individuals based on the known releases from a year/zone
      {
        if(release_sea(s2-1)!=1&&release_sea(s2-1)!=2&&release_sea(s2-1)!=3&&release_sea(s2-1)!=4&&release_sea(s2-1)!=11&&release_sea(s2-1)!=12)//For the non-winter months apply the survival AND movement
        {
          Ntagmat(s,z,s2)=migmat(release_sea(s2-1)-4)*elem_prod(exp(-Z(s2-1)),Ntagmat(s,z,s2-1)); //matrix multiplication of previous season's N, survival, and movement
          //cout<<"Migmat:"<<migmat(release_sea(s2-1)-2)<<endl<<endl;
        }
        if(release_sea(s2-1)==1||release_sea(s2-1)==2||release_sea(s2-1)==3||release_sea(s2-1)==4||release_sea(s2-1)==11||release_sea(s2-1)==12)//For the non-winter months apply the survival only
        {
          Ntagmat(s,z,s2)=elem_prod(exp(-Z(s2-1)),Ntagmat(s,z,s2-1)); //matrix multiplication of previous season's N and survival
        }
      }
      
    }
  }
  
  //cout<<"Number Tagged Matrix:"<<endl<<Ntagmat<<endl;
  //exit(1);

FUNCTION calc_negLL

  negLLrecapture=0.0; //Set the initial likelihood values to 0
  negLLq=0.0;
  negLLmig=0.0;
  negLL=0.0;
  
  for(s=1;s<=numsts;s++)
  {
    for(z=1;z<=numzone;z++)
    {
      likeli(s,z)=0.0; //setting the initial likelihood to 0
      expected(s,z)=0.0; //setting the initial number of tagged invididuals in subsequent years to 0
      
      if(releases(s,z)>0)//if there are actual releases from a certain month/region AND if there is effort in the recapture region, the likelihood gets calculated for the first time step
      {
        if(Ntagmat(s,z,s,z)*pobs(s,z)>0)
        {
          expected(s,z,s,z)=Ntagmat(s,z,s,z)*pobs(s,z); //calculated the expected returns, from the number of individuals and the probability of observation

          //likelihood function NEGATIVE BINOMIAL:
          likeli(s,z,s,z)=gammln(k)+gammln(recapture(s,z,s,z)+1.0)+k*log(k+expected(s,z,s,z))+recapture(s,z,s,z)*log(k+expected(s,z,s,z)); //likelihood is calculated on two lines because of programming issues
          likeli(s,z,s,z)=likeli(s,z,s,z)-(gammln(k+recapture(s,z,s,z))+k*log(k)+recapture(s,z,s,z)*log(expected(s,z,s,z)));
          
        }
        for(s2=s+1;s2<=numsts;s2++) //cycles through all subsequent time steps and regions and calculates the likelihood 
        {
          for(z2=1;z2<=numzone;z2++)
          {
            if(Ntagmat(s,z,s2,z2)*pobs(s2,z2)>0)
            {
              expected(s,z,s2,z2)=Ntagmat(s,z,s2,z2)*pobs(s2,z2);  //calculated the expected returns, from the number of individuals and the probability of observation

              //likelihood function NEGATIVE BINOMIAL:
              likeli(s,z,s2,z2)=gammln(k)+gammln(recapture(s,z,s2,z2)+1.0)+k*log(k+expected(s,z,s2,z2))+recapture(s,z,s2,z2)*log(k+expected(s,z,s2,z2)); //likelihood is calculated on two lines because of programming issues
              likeli(s,z,s2,z2)=likeli(s,z,s2,z2)-(gammln(k+recapture(s,z,s2,z2))+k*log(k)+recapture(s,z,s2,z2)*log(expected(s,z,s2,z2)));

              //negLL+=likeli(s,z,s2,z2);
              //cout<<"Position: "<<s<<" "<<z<<" "<<s2<<" "<<z2<<" ";
              //cout<<"Ntagmat: "<<Ntagmat(s,z,s2,z2) << " ";
              //cout<<"Recapture Number: "<<recapture(s,z,s2,z2)<<" ";
              //cout<<"Likelihood: "<<likeli(s,z,s2,z2)<<endl;
            }
          }
        }
        negLLrecapture+=sum(likeli(s,z)); //sums the likelihood of each individual observation into a single value
      }
      //cout<<"Neg Log Likelihood: "<<negLL<<endl;
    }
  }
  //exit(1);

  for(s=1;s<=numsea-7;s++) //for each of the MONTHLY movement matrix values, cycle through and sum the effect of the Dirichlet distribution prior
  {
    negLLmig+=-sum(elem_prod(effsamplesize*migrationprior(1),log(migmat(s)))); 
  }
  negLLmig+=-sum(elem_prod(effsamplesize*migrationprior(2),log(migmat(6)))); //for the HALF-YEARLY movement matrix values, cycle through and sum the effect of the Dirichlet distribution prior

  //calculates the average log adjusted total catchability based on just the non-zero values
  sum_qtot=0;
  count_qtot=0;
  for(s=1;s<=numsea;s++)
  {
    for(z=1;z<=numzone;z++)
    {
      if(qtot(s,z)!=0.0) //if there is an actual q value, count it and sum the values
      {
      sum_qtot+=log(qtot(s,z));
      //sum_qtot+=qtot(s,z);
      count_qtot+=1.0;
      }
    }
  }
  mean_qtot=sum_qtot/count_qtot; //calculate the average catchability based on the sum and the count

  //cout<<"Mean Total Catchability Value: "<<mean_qtot<<endl<<endl;
  //cout<<"Catchability: "<<qtot<<endl;
  //cout<<"Likelihood before q: "<<negLL<<endl<<endl;

  //calculates the deviation from the mean of catchability values and adds them to prior penalty component
  for(s=1;s<=numsea;s++)
  {
    for(z=1;z<=numzone;z++)
    {
      if(qtot(s,z)!=0.0)
      {
      negLLq+=0.5*square(log(qtot(s,z))-mean_qtot)/square(1.0); //calculates the deviation from the mean catchability value, the value in the denominator is the variance of assumed distribution
      }
    }
  }
  negLL=negLLrecapture+negLLmig+negLLq; //sum likelihood and penalty components into the objective function

  //cout<<"Likelihood after q: "<<negLL<<endl<<endl;
  //cout<<"qtot:"<<qtot<<endl<<endl;
  //exit(1);

FUNCTION calc_mageff
  est_mageff(1,1)=0;
  est_mageff(2,1)=0;
  est_mageff(3,1)=0;
  est_mageff(4,1)=0;
  est_mageff(5,1)=0;
  est_mageff(6,1)=0;
  est_mageff(7,1)=0;
  est_mageff(8,1)=0;
  est_mageff(9,1)=0;
  est_mageff(10,1)=0;
  est_mageff(11,1)=0;
  est_mageff(12,1)=exp(log_mageffv(1));
  est_mageff(13,1)=exp(log_mageffv(2));
  est_mageff(14,1)=exp(log_mageffv(3));
  est_mageff(15,1)=exp(log_mageffv(4));
  est_mageff(16,1)=exp(log_mageffv(5));
  est_mageff(17,1)=0;
  est_mageff(18,1)=0;
  est_mageff(19,1)=0;
  est_mageff(20,1)=0;
  est_mageff(21,1)=0;
  est_mageff(22,1)=0;
  est_mageff(23,1)=mageff_bound(1);  //exp(log_mageffv(6));
  est_mageff(24,1)=exp(log_mageffv(6));
  est_mageff(25,1)=exp(log_mageffv(7));
  est_mageff(26,1)=exp(log_mageffv(8));
  est_mageff(27,1)=exp(log_mageffv(9));
  est_mageff(28,1)=exp(log_mageffv(10));
  est_mageff(29,1)=0;
  est_mageff(30,1)=0;
  est_mageff(31,1)=0;
  est_mageff(32,1)=0;
  est_mageff(33,1)=0;
  est_mageff(34,1)=0;
  est_mageff(35,1)=0;
  est_mageff(36,1)=exp(log_mageffv(11));
  est_mageff(37,1)=exp(log_mageffv(12));
  est_mageff(38,1)=exp(log_mageffv(13));
  est_mageff(39,1)=exp(log_mageffv(14));
  est_mageff(40,1)=exp(log_mageffv(15));
  est_mageff(41,1)=0;
  est_mageff(42,1)=0;


  est_mageff(1,2)=0;
  est_mageff(2,2)=0;
  est_mageff(3,2)=0;
  est_mageff(4,2)=0;
  est_mageff(5,2)=0;
  est_mageff(6,2)=0;
  est_mageff(7,2)=0;
  est_mageff(8,2)=0;
  est_mageff(9,2)=0;
  est_mageff(10,2)=0;
  est_mageff(11,2)=exp(log_mageffv(16));
  est_mageff(12,2)=exp(log_mageffv(17));
  est_mageff(13,2)=exp(log_mageffv(18));
  est_mageff(14,2)=exp(log_mageffv(19));
  est_mageff(15,2)=exp(log_mageffv(20));
  est_mageff(16,2)=exp(log_mageffv(21));
  est_mageff(17,2)=exp(log_mageffv(22));
  est_mageff(18,2)=0;
  est_mageff(19,2)=0;
  est_mageff(20,2)=0;
  est_mageff(21,2)=0;
  est_mageff(22,2)=exp(log_mageffv(23));
  est_mageff(23,2)=exp(log_mageffv(24));
  est_mageff(24,2)=exp(log_mageffv(25));
  est_mageff(25,2)=exp(log_mageffv(26));
  est_mageff(26,2)=exp(log_mageffv(27));
  est_mageff(27,2)=exp(log_mageffv(28));
  est_mageff(28,2)=exp(log_mageffv(29));
  est_mageff(29,2)=exp(log_mageffv(30));
  est_mageff(30,2)=0;
  est_mageff(31,2)=0;
  est_mageff(32,2)=0;
  est_mageff(33,2)=0;
  est_mageff(34,2)=exp(log_mageffv(31));
  est_mageff(35,2)=exp(log_mageffv(32));
  est_mageff(36,2)=exp(log_mageffv(33));
  est_mageff(37,2)=exp(log_mageffv(34));
  est_mageff(38,2)=exp(log_mageffv(35));
  est_mageff(39,2)=exp(log_mageffv(36));
  est_mageff(40,2)=exp(log_mageffv(37));
  est_mageff(41,2)=exp(log_mageffv(38));
  est_mageff(42,2)=0;

  est_mageff(1,3)=exp(log_mageffv(39));
  est_mageff(2,3)=exp(log_mageffv(40));
  est_mageff(3,3)=exp(log_mageffv(41));
  est_mageff(4,3)=exp(log_mageffv(42));
  est_mageff(5,3)=mageff_bound(2);  //exp(log_mageffv(43));
  est_mageff(6,3)=exp(log_mageffv(43));
  est_mageff(7,3)=exp(log_mageffv(44));
  est_mageff(8,3)=0;
  est_mageff(9,3)=exp(log_mageffv(45));
  est_mageff(10,3)=0;
  est_mageff(11,3)=exp(log_mageffv(46));
  est_mageff(12,3)=exp(log_mageffv(47));
  est_mageff(13,3)=exp(log_mageffv(48));
  est_mageff(14,3)=exp(log_mageffv(49));
  est_mageff(15,3)=exp(log_mageffv(50));
  est_mageff(16,3)=exp(log_mageffv(51));
  est_mageff(17,3)=exp(log_mageffv(52));
  est_mageff(18,3)=mageff_bound(2); //exp(log_mageffv(50));
  est_mageff(19,3)=0;
  est_mageff(20,3)=0;
  est_mageff(21,3)=exp(log_mageffv(53));
  est_mageff(22,3)=exp(log_mageffv(54));
  est_mageff(23,3)=exp(log_mageffv(55));
  est_mageff(24,3)=exp(log_mageffv(56));
  est_mageff(25,3)=exp(log_mageffv(57));
  est_mageff(26,3)=exp(log_mageffv(58));
  est_mageff(27,3)=exp(log_mageffv(59));
  est_mageff(28,3)=mageff_bound(2);     //exp(log_mageffv(63));
  est_mageff(29,3)=mageff_bound(2);     //exp(log_mageffv(64));
  est_mageff(30,3)=mageff_bound(2);   //exp(log_mageffv(57));
  est_mageff(31,3)=mageff_bound(1);    // exp(log_mageffv(58));
  est_mageff(32,3)=exp(log_mageffv(60));
  est_mageff(33,3)=exp(log_mageffv(61));
  est_mageff(34,3)=0;
  est_mageff(35,3)=exp(log_mageffv(62));
  est_mageff(36,3)=exp(log_mageffv(63));
  est_mageff(37,3)=exp(log_mageffv(64));
  est_mageff(38,3)=exp(log_mageffv(65));
  est_mageff(39,3)=exp(log_mageffv(66));
  est_mageff(40,3)=exp(log_mageffv(67));
  est_mageff(41,3)=exp(log_mageffv(68));
  est_mageff(42,3)=exp(log_mageffv(69));


  est_mageff(1,4)=0;
  est_mageff(2,4)=0;
  est_mageff(3,4)=0;
  est_mageff(4,4)=0;
  est_mageff(5,4)=0;
  est_mageff(6,4)=0;
  est_mageff(7,4)=0;
  est_mageff(8,4)=0;
  est_mageff(9,4)=0;
  est_mageff(10,4)=mageff_bound(2);    //exp(log_mageffv(69));
  est_mageff(11,4)=mageff_bound(2);   //exp(log_mageffv(70));
  est_mageff(12,4)=exp(log_mageffv(70));
  est_mageff(13,4)=exp(log_mageffv(71));
  est_mageff(14,4)=exp(log_mageffv(72));
  est_mageff(15,4)=mageff_bound(1);   //exp(log_mageffv(74));
  est_mageff(16,4)=exp(log_mageffv(73));
  est_mageff(17,4)=mageff_bound(2);  //exp(log_mageffv(76));
  est_mageff(18,4)=0;
  est_mageff(19,4)=0;
  est_mageff(20,4)=0;
  est_mageff(21,4)=0;
  est_mageff(22,4)=exp(log_mageffv(74));
  est_mageff(23,4)=exp(log_mageffv(75));
  est_mageff(24,4)=exp(log_mageffv(76));
  est_mageff(25,4)=exp(log_mageffv(77));
  est_mageff(26,4)=exp(log_mageffv(78));
  est_mageff(27,4)=exp(log_mageffv(79));
  est_mageff(28,4)=exp(log_mageffv(80));
  est_mageff(29,4)=mageff_bound(2);    //exp(log_mageffv(92));
  est_mageff(30,4)=0;
  est_mageff(31,4)=0;
  est_mageff(32,4)=0;
  est_mageff(33,4)=0;
  est_mageff(34,4)=exp(log_mageffv(81));
  est_mageff(35,4)=exp(log_mageffv(82));
  est_mageff(36,4)=exp(log_mageffv(83));
  est_mageff(37,4)=exp(log_mageffv(84));
  est_mageff(38,4)=mageff_bound(1);     //exp(log_mageffv(97));
  est_mageff(39,4)=mageff_bound(1);     //exp(log_mageffv(98));
  est_mageff(40,4)=mageff_bound(1);     //exp(log_mageffv(99));
  est_mageff(41,4)=exp(log_mageffv(85));
  est_mageff(42,4)=0;
  
 

FUNCTION set_solutions
  migmatsolutions1=migmat(1); //set the movement matrices to SD report variables
  migmatsolutions2=migmat(2);
  migmatsolutions3=migmat(3);
  migmatsolutions4=migmat(4);
  migmatsolutions5=migmat(5);
  migmatsolutions6=migmat(6);
  
  qtotsolutions(1)=qtot(5,1); //set the catchability values to SD report variables
  qtotsolutions(2)=qtot(6,1);
  qtotsolutions(3)=qtot(7,1);
  qtotsolutions(4)=qtot(8,1);
  qtotsolutions(5)=qtot(9,1);
  qtotsolutions(6)=qtot(10,1);
  qtotsolutions(7)=qtot(11,1);
  qtotsolutions(8)=qtot(5,2);
  qtotsolutions(9)=qtot(6,2);
  qtotsolutions(10)=qtot(7,2);
  qtotsolutions(11)=qtot(8,2);
  qtotsolutions(12)=qtot(9,2);
  qtotsolutions(13)=qtot(10,2);
  qtotsolutions(14)=qtot(11,2);
  qtotsolutions(15)=qtot(12,2);
  qtotsolutions(16)=qtot(1,3);
  qtotsolutions(17)=qtot(4,3);
  qtotsolutions(18)=qtot(5,3);
  qtotsolutions(19)=qtot(6,3);
  qtotsolutions(20)=qtot(7,3);
  qtotsolutions(21)=qtot(8,3);
  qtotsolutions(22)=qtot(9,3);
  qtotsolutions(23)=qtot(10,3);
  qtotsolutions(24)=qtot(11,3);
  qtotsolutions(25)=qtot(12,3);
  qtotsolutions(26)=qtot(4,4);
  qtotsolutions(27)=qtot(5,4);
  qtotsolutions(28)=qtot(6,4);
  qtotsolutions(29)=qtot(7,4);
  qtotsolutions(30)=qtot(8,4);
  qtotsolutions(31)=qtot(9,4);
  qtotsolutions(32)=qtot(10,4);
  qtotsolutions(33)=qtot(11,4);


FUNCTION mcmcoutput
  //this section outputs the results of the MCMC output into a single line of values
  ofstream mcmcout("Coston_ME85params.mcmc",ios::app);
  {
    mcmcout<<negLL<<" "<<log_M<<log_q<<theta<<" "<<log_k;
    
    for(i=1;i<=numsea-6;i++)
    {
      for(z=1;z<=numzone-1;z++)
      {
        mcmcout<<logmigparam(i,z);
      }
    }

    mcmcout<<" "<<M<<" "<<Myr<<" "<<k<<" ";
    for(i=1;i<=numzone+numtheta;i++)
    {
      mcmcout<<qtotsolutions(i)<<" ";
    }
    
    for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        mcmcout<<F(s,z)<<" ";
      }
    }
    
    for(i=1;i<=numsea-6;i++)
    {
      for(z=1;z<=numzone;z++)
      {
        mcmcout<<migmat(i,z)<<" ";
      }
    }


 // only output the non zero log_mageffv
     for (i=1;i<=nummageff;i++)
     {
     mcmcout<<log_mageffv(i)<<" ";
     }
  
    mcmcout<<endl;    
  }

FUNCTION write_waic_row

  // Writes one row of pointwise log-likelihoods (NOT negative) for WAIC.
  // Column order: for s=1..numsts, z=1..numzone, s2=1..numsts, z2=1..numzone.
  // Header is created automatically on first run.

  // Decide whether to write header
  bool write_header = false;
  {
    ifstream chk("Coston_waic_loglik.csv");
    if (!chk) write_header = true;
  }

  ofstream w(write_header ? "Coston_waic_loglik.csv" : "Coston_waic_loglik.csv", 
             write_header ? ios::out : ios::app);

  w.setf(std::ios::fixed);
  w.precision(10);

  // Header with column labels (once)
  if (write_header)
  {
    bool first = true;
    for (int s_rel=1; s_rel<=numsts; ++s_rel)
      for (int z_rel=1; z_rel<=numzone; ++z_rel)
        for (int s_rec=1; s_rec<=numsts; ++s_rec)
          for (int z_rec=1; z_rec<=numzone; ++z_rec)
          {
            if (first) { first = false; }
            else       { w << ","; }
            w << "s"  << s_rel
              << "_z" << z_rel
              << "_s2"<< s_rec
              << "_z2"<< z_rec;
          }
    w << "\n";
  }

  // Data row: log-likelihood = -likeli(...)
  {
    bool first = true;
    for (int s_rel=1; s_rel<=numsts; ++s_rel)
      for (int z_rel=1; z_rel<=numzone; ++z_rel)
        for (int s_rec=1; s_rec<=numsts; ++s_rec)
          for (int z_rec=1; z_rec<=numzone; ++z_rec)
          {
            // Note: likeli holds NEGATIVE log-likelihood; WAIC needs log-likelihood
            double ll = -value(likeli(s_rel, z_rel, s_rec, z_rec));
            if (first) { w << ll; first = false; }
            else       { w << "," << ll; }
          }
    w << "\n";
  }

  w.close();



REPORT_SECTION

     cout<<"expected  Matrix:"<<endl<<sum(expected)<<endl;
      cout<<"releases  Matrix:"<<endl<<sum(releases)<<endl;
       cout<<"recapture  Matrix:"<<endl<<sum(recapture)<<endl;
       cout<<"Number Tag Out Matrix:"<<endl<<sum(Ntagmat)<<endl;
 for(s=1;s<=42;s++)
  {
    for(z=1;z<=numzone;z++)
    {
      summ=0.0;
      summ0=0.0;
     // summ0=0.0;
       for(s2=1;s2<=42;s2++)
         {
         for (z2=1;z2<=numzone;z2++)
          {
         summ+=expected(s2,z2,s,z);
         summ0+=recapture(s2,z2,s,z);
        // summ0+=expected0(s2,z2,s,z);
          }
         }
       catches(s,z)=summ;
       orecapture(s,z)=summ0;
      // catch0(s,z)=summ0;
       
     }
  }

// calculate the ratios of recapture and expected0
 for(s=1;s<=42;s++)
  {
    for(z=1;z<=numzone;z++)
    {
      summ=0.0;
       for(s2=1;s2<=42;s2++)
         {
         for (z2=1;z2<=numzone;z2++)
          {
         summ+=Ntagmat(s2,z2,s,z);
          }
         }
       N_pop(s,z)=summ;
     }
  }

    //prints out the results of a single run in the same format at the MCMC output

    //report<<negLL<<" "<<log_M<<log_q<<theta<<" "<<log_k;
    report<<negLL<<" "<<M<<" "<<Myr;
    report<<" "<<sum(expected)<<" ";
    report<<sum(releases)<<" ";
    report<<sum(recapture)<<" ";
    report<<sum(N_pop)<<" ";

    for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        report<<F(s,z)<<" ";
      }
    }

   for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        report<<effort(s,z)<<" ";
      }
    }
    
   for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        report<<Z(s,z)<<" ";
      }
    }
  for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        report<<N_pop(s,z)<<" ";
      }
    }
  for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        report<<catches(s,z)<<" ";
      }
    }

  for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        report<<orecapture(s,z)<<" ";
      }
    }

  for(z=1;z<=numzone;z++)
    {
      for(s=1;s<=numsts;s++)
      {
        report<<est_mageff(s,z)<<" ";
      }
    }


    for(i=1;i<=numsea-6;i++)
    {
      for(z=1;z<=numzone-1;z++)
      {
        report<<logmigparam(i,z)<<" ";
      }
    }

    
    for(i=1;i<=numzone+numtheta;i++)
    {
      report<<qtotsolutions(i)<<" ";
    }
    
    
    for(i=1;i<=numsea-6;i++)
    {
      for(z=1;z<=numzone;z++)
      {
        report<<migmat(i,z)<<" ";
      }
    }

   report<<log_M<<" "<<log_q<<theta<<" "<<log_k<<" "<<k<<" ";
   report<<nummageff<<" ";

    report<<endl;


