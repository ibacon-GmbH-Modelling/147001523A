/*
  FILE: test_derivatives.cpp version of 20220128
  for BYOM_v6/DEBtox2019_v45b
 
 Below: all licences and copyright notices of the code used here.
 
======================
 
 Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

 Copyright 2010-2012 Karsten Ahnert
 Copyright 2011-2013 Mario Mulansky
 Copyright 2013 Pascal Germroth
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 
 =====================

 % * Author: Tjalling Jager
 % * Date: September 2020
 % * Web support: <http://www.debtox.info/byom.html>
 % * Back to index <walkthrough_debtox2019.html>

 %  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
 %  This source code is licensed under the MIT-style license found in the
 %  LICENSE.txt file in the root directory of BYOM. 
 
 ======================

 Edits to apply it to the problem at hand by Dr. Carlo Romoli - ibacon GmbH

 Some of the technical solutions in the code have been taken from the example 
 reported in the Boost-libraries documentation to solve differential
 equations (see above copyright notice).

 This C++ code is a translation of the DEBtox2019 MATLAB code 
 developed by Dr. Tjalling Jager (see above copyright notice). The equations
 are those reported in the derivatives.m and read_scen.m files of the 
 DEBtox2019 package.

 The connection between C++ and MATLAB has been done using the 
 MATLAB C++ MEX APIs
 
 =======================
 */


#include <iostream>
#include <vector>
#include <algorithm>

#include <boost/numeric/odeint.hpp>

#include "mex.hpp"
#include "mexAdapter.hpp"

using matlab::mex::ArgumentList;
using namespace matlab::data;
using namespace matlab::mex;

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

// call a global pointer to start the matlab engine so that it is always
// available and does not need to be called all the time
//std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr = connectMATLAB();

class DEBderi {
	// the parameters were originally in a structure.
	// this has been converted into vectors for easieness and performance
	std::vector<std::vector<double>> vectors;  // parameters of the model
	std::vector<double> scalars;               // more global parameters of the model
    double ci;                                 // external concentration that apparently is a double
    matlab::data::TypedArray<double> int_coll; // 2D vector containing the scenario table
	std::vector<double> int_coll_times;        // time vector of the scenario table
    int int_type;                              // type of scenario (2 for constant, 4 for linear interpolation)
    std::vector<double> timevar;               // 2 element array telling if we have a variable profile or not
    std::shared_ptr<matlab::engine::MATLABEngine> mateng; // MATLAB engine for debugging purposes (to allow printouts)
    
	public:
		DEBderi(std::vector<double> scalar_pars,
                std::vector<std::vector<double>> vector_pars,
                double conc,
                matlab::data::TypedArray<double> int_coll_arr,
				std::vector<double> int_coll_times_arr,
                int int_type_val,
                std::vector<double> timevar_arr,
                std::shared_ptr<matlab::engine::MATLABEngine> mateng2) : scalars(scalar_pars), // initializer list
                                                                         vectors(vector_pars),
                                                                         ci(conc),
                                                                         int_coll(int_coll_arr),
							                                             int_coll_times(int_coll_times_arr),
                                                                         int_type(int_type_val),
                                                                         timevar(timevar_arr),
                                                                         mateng(mateng2){}
		void operator() ( state_type &x , state_type &dxdt , const double t ) // not declaring x as constant otherwise bad?
		{
			/* insert all the derivatives from the DEB model */
            // unpack parameters. They will need to be passed in the same
            // order by the main function. Unless we generate a new structure
            // we will not have the name of the fields available

            // The parameters are passed through a C++ vector to this class
            // in order to increase speed. Reading the MATLAB object in the
            // class with the derivatives would be to heavy
            double FBV = scalars[0];
            double KRV = scalars[1];     // part. coeff. repro buffer and structure (kg/kg)
            double kap = scalars[2];     // approximation for kappa (-)
            double yP  = scalars[3];      // product of yVA and yAV (-)

            double L0   = scalars[4];   // body length at start (mm)
            double Lp   = scalars[5];   // body length at puberty (mm)
            double Lm   = scalars[6];   // maximum body length (mm)
            double rB   = scalars[7];   // von Bertalanffy growth rate constant (1/d)
            double Rm   = scalars[8];   // maximum reproduction rate (#/d)
            double f    = scalars[9];    // scaled functional response (-)
            double hb   = scalars[10];   // background hazard rate (d-1)
 
            // unpack extra parameters for specific cases
            double Lf   = scalars[11];   // body length at half-saturation feeding (mm)
            double Tlag = scalars[12]; // lag time for start development (d)
            // unpack model parameters for the response to toxicants
            double kd   = scalars[13];   // dominant rate constant (d-1)
            double zb   = scalars[14];   // effect threshold energy budget ([C])
            double bb   = scalars[15];   // effect strength energy-budget effects (1/[C])
            double zs   = scalars[16];   // effect threshold survival ([C])
            double bs   = scalars[17];   // effect strength survival (1/([C] d))

            double Lj = scalars[18]; // length at metamorphosis (for abj models) No need for Daphnia
            double Lm_ref = scalars[19];
			double MF = scalars[20];
			double a = scalars[21];
			
			hb = a * std::pow(hb,a) * std::pow(t,(a-1)); // option for Weibull mortalty when a is not 1

            std::vector<double> feedbacks(vectors[0]);
            std::vector<double> moa(vectors[1]);

            // initial conditions read from the input and set so that they
            // do not become negative (as in original code)
            x[0]=std::max(x[0],0.);
            x[1]=std::max(x[1],0.);
            x[2]=std::max(x[2],0.);
            x[3]=std::max(x[3],0.);
            
            double c=ci;  // concentration or concentration scenario
            
            // only in case we have variable concentrations
			if ((int)timevar[0] == 1){
                //stream << "calling external function\n";
                //displayOnMATLAB(stream);
                c = read_scen(ci, t, MF, int_coll, int_coll_times, int_type, timevar); // for time varying concentrations
            }

            x[1] = std::max(1e-3 * L0, x[1]);
            
            if (Lf > 0){
                f = f / (1 + (Lf * Lf * Lf)/(x[1] * x[1] * x[1])); // hyperbolic relationship for f with body volume
            }
            if (Lj > 0) {// to include acceleration until metamorphosis ...
                f = f * std::min(1.,x[1]/Lj); // this implies lower f for L<Lj
            }
            
            double s = bb*std::max(0.,x[0]-zb); // stress level for metabolic effects
            double h = bs*std::max(0.,x[0]-zs); // hazard rate for effects on survival
			
			h = std::min(111.,h);  // maximise the hazard rate to 99% mortality in 1 hour
			// Note: this helps in extreme conditions, as the system becomes stiff for
            // very high hazard rates. This is especially needed for EPx calculations,
            // where the MF is increased until there is effect on all endpoints!
            
            // 5 MODE OF ACTION
            double sA = std::min(1.,moa[0] * s); // assimilation/feeding (maximise to 1 to avoid negative values for 1-sA)
            double sM = moa[1] * s;              // maintenance (somatic and maturity)
            double sG = moa[2] * s;              // growth costs
            double sR = moa[3] * s;              // reproduction costs
            double sH = moa[4] * s;              // also include hazard to reproduction

            dxdt[1] = rB * ((1+sM)/(1+sG)) * (f*Lm*((1-sA)/(1+sM)) - x[1]); // ODE for body length
            
            double fR = f; // if there is no starvation, f for reproduction is the standard f
            // starvation rules can modify the outputs here
            if (dxdt[1] < 0){ // then we are looking at starvation and need to correct things
                fR = (f - kap * (x[1]/Lm) * ((1+sM)/(1-sA)))/(1-kap); // new f for reproduction alone
                if (fR >= 0){  // then we are in the first stage of starvation: 1-kappa branch can help pay maintenance
                    dxdt[1] = 0; // stop growth, but don't shrink
                } else {        // we are in stage 2 of starvation and need to shrink to pay maintenance
                    fR = 0; // nothing left for reproduction
                    dxdt[1] = (rB*(1+sM)/yP) * ((f*Lm/kap)*((1-sA)/(1+sM)) - x[1]); // shrinking rate
                }
            }

            double R  = 0; // reproduction rate is zero, unless ... 
            if (x[1] >= Lp){ // if we are above the length at puberty, reproduce
                //R = std::max(0.,(Rm/(1+sR)) * (fR*Lm*(x[1]*x[1])*(1-sA) - (Lp*Lp*Lp)*(1+sM))/(Lm*Lm*Lm - Lp*Lp*Lp));
                R = std::max(0.,(exp(-sH)*Rm/(1+sR)) * (fR*Lm*(x[1]*x[1])*(1-sA) - (Lp*Lp*Lp)*(1+sM))/(Lm*Lm*Lm - Lp*Lp*Lp));
            }
            dxdt[2] = R;                 // cumulative reproduction rate
            dxdt[3]  = -(h + hb) * x[3]; // change in survival probability (incl. background mort.)

            // For the damage dynamics, there are four feedback factors x* that obtain a
            // value based on the settings in the configuration vector glo.feedb: a
            // vector with switches for various feedbacks: [surface:volume on uptake,
            // surface:volume on elimination, growth dilution, losses with
            // reproduction].

            // this operation has to be handled with care
            // element-wise product
            feedbacks[0] = feedbacks[0] * Lm_ref/x[1];
            feedbacks[1] = feedbacks[1] * Lm_ref/x[1];
            feedbacks[2] = feedbacks[2] * (3/x[1])*dxdt[1];
            feedbacks[3] = feedbacks[3] * R*FBV*KRV;
            
            //double xu = std::max(1.,feedbacks[0]); // if switch for surf:vol scaling is zero, the factor must be 1 and not 0!
            //double xe = std::max(1.,feedbacks[1]); // if switch for surf:vol scaling is zero, the factor must be 1 and not 0!
            double xu = feedbacks[0];
            if (feedbacks[0] == 0) {xu = 1;}
            double xe = feedbacks[1];
            if (feedbacks[1] == 0) {xe = 1;}
            double xG = feedbacks[2];              // factor for growth dilution
            double xR = feedbacks[3];              // factor for losses with repro

            xG = std::max(0.,xG); 
            // NOTE NOTE: reverse growth dilution (concentration by shrinking) is now
            // turned OFF as it leads to runaway situations that lead to failure of the
            // ODE solvers. However, this needs some further thought!
            dxdt[0] = kd * (xu * c - xe * x[0]) - (xG + xR) * x[0]; // ODE for scaled damage

            if (x[1] <= 0.5 * L0){ // if an animal has size less than half the start size ...
                dxdt[1] = 0.; // don't let it grow or shrink any further (to avoid numerical issues)
            }
            
            if (t<Tlag){
                //derivatives are non-zero only if time is greater than Tlag
                dxdt[0] = 0;
                dxdt[1] = 0;
                dxdt[2] = 0;
                dxdt[3] = 0;
            }
	    }

        double read_scen(double c, double t, double MF, matlab::data::TypedArray<double> int_coll,
		                 std::vector<double> int_coll_times, int int_type, std::vector<double> timevar){
            // function copied from DEBtox to avoid calling matlab code from here
            // it copied only the case with -1, the case called from the derivatives file
			
            double out_c=0;
            int size_t_int_coll = (int)int_coll.getDimensions()[0];
            switch (int_type){
                case 1:
                    {
                    // should not be needed.
                        break;
                    }
                case 2:
                    {
                    if (timevar.size()==2 && timevar[1] > 0){
                        out_c = MF * int_coll[(int)timevar[1]-1][1];
                    }
                    else{
                        // substitute the find function here with a C++ implementation of this case
                        // should not be needed now, but who knows
                        //out_c = MF * int_coll(find(int_coll(:,1)<=t,1,'last'),2);
                        // use std::find_if increases the performance massively!
                        auto it = std::find_if(int_coll_times.rbegin(),int_coll_times.rend(),[&](const double& i){return i<=t;});
						int ii = it - int_coll_times.rbegin();
                        ii = int_coll_times.size() - 1 - ii;
                        out_c = MF * int_coll[ii][1];    
                    }
                    break;
                    }
                case 3:
                    {
                    double kc = int_coll[size_t_int_coll-1][1];
                    if (timevar.size()==2 && timevar[1] > 0){
                        int ind_i = (int)timevar[1]-1;
                        double c0 = int_coll[ind_i][1];
                        double t0 = t - int_coll[ind_i][0];
                        out_c = MF * c0 * exp(-kc*t0);
                    }
                    else{
                        auto it = std::find_if(int_coll_times.rbegin(),int_coll_times.rend(),[&](const double& i){return i<=t;});
						int ii = it - int_coll_times.rbegin();
                        ii = int_coll_times.size() - 1 - ii;
                        double c0 = int_coll[ii][1];
                        double t0 = t - int_coll[ii][0];
                        out_c = MF * c0 * exp(-kc*t0);                    
                        }
                    break;
                    }
                case 4:
                    {
                    if (timevar.size()==2 && timevar[1] > 0){
                        int ind_i = (int)timevar[1]-1;   
                        out_c = int_coll[ind_i][1] * MF + (t - int_coll[ind_i][0]) * int_coll[ind_i][2] * MF;
                    }
                    else{
						auto it = std::find_if(int_coll_times.rbegin(),int_coll_times.rend(),[&](const double& i){return i<=t;});
						int ii = it - int_coll_times.rbegin();
                        ii = int_coll_times.size() - 1 - ii;
                        out_c = int_coll[ii][1] * MF + (t-int_coll[ii][0]) * int_coll[ii][2] * MF;
                    }
                    break;
                    }
            }
            return out_c;
        }

        void displayOnMATLAB(std::ostringstream& stream) {
			// function to printout stuff.
			// Work on a wat to make the inheritance instead of
			// duplicating code
            ArrayFactory factory;
            // Pass stream content to MATLAB fprintf function
            mateng->feval(u"fprintf", 0,
              std::vector<Array>({ factory.createScalar(stream.str()) }));
            // Clear stream buffer
            stream.str("");
        }
};

//[ integrate_observer
// structure containing the function to store the states at each step
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
//]

class MexFunction : public matlab::mex::Function { 
    // create pointer to matlab engine
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr2 = getEngine();
    // Factory to create MATLAB data arrays
    ArrayFactory factory;
    public:
      // Print strings during exectution. Useful for DEBUG
      void displayOnMATLAB(std::ostringstream& stream) {
          // Pass stream content to MATLAB fprintf function
          matlabPtr2->feval(u"fprintf", 0,
              std::vector<Array>({ factory.createScalar(stream.str()) }));
          // Clear stream buffer
          stream.str("");
      }

      void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs){    
          using namespace std;
          using namespace boost::numeric::odeint;
          // Create an output stream
          // ostringstream stream; // needed in case of needing DEBUG

          /* input parameters in the same order as they appear in derivatives.m
           * -time range, 
           * -initial conditions
           * -par
           * -c
           * -glo
           * -initial step (intial dt for the solver)
		   * -abstol (error tolerances of the ODE solver)
		   * -reltol
		   * -max step size
		   * (the maximum step size is needed to avoid the dense adaptive stepper to 
		   * perform steps that are too large)
           */

          // time range
          matlab::data::TypedArray<double> inArray = inputs[0];
          vector<double> time_vector(inArray.begin(), inArray.end());
          // initial conditions
          matlab::data::TypedArray<double> inArray2 = inputs[1];
          vector<double> init_states(inArray2.begin(), inArray2.end());
          // par structure
          matlab::data::StructArray inStructArrayPar = inputs[2];
          // concentration c
          double conc = inputs[3][0];
          // glo structure
	      matlab::data::StructArray inStructArrayGlo = inputs[4];
	      // additional argument as I need an initial dt to be passed for the solver
          double dt = inputs[5][0];
          double AbsErr = inputs[6][0]; // tolerances for the ODE solver
          double RelErr = inputs[7][0];
		  double MaxStep = inputs[8][0]; // maximum step-size

          vector<double> scalar_pars;
          vector<std::vector<double>> vector_pars;

          // Extract all the parameters from glo and par
          matlab::data::Array tempconv = inStructArrayGlo[0]["FBV"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayGlo[0]["KRV"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayGlo[0]["kap"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayGlo[0]["yP"];
          scalar_pars.push_back(tempconv[0]);

          // unpack model parameters for the basic life history
          tempconv = inStructArrayPar[0]["L0"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["Lp"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["Lm"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["rB"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["Rm"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["f"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["hb"];
          scalar_pars.push_back(tempconv[0]);

          // unpack extra parameters for specific cases
          tempconv = inStructArrayPar[0]["Lf"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["Tlag"];
          scalar_pars.push_back(tempconv[0]);

          // unpack model parameters for the response to toxicants
          tempconv = inStructArrayPar[0]["kd"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["zb"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["bb"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["zs"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayPar[0]["bs"];
          scalar_pars.push_back(tempconv[0]);

          tempconv = inStructArrayPar[0]["Lj"];
          scalar_pars.push_back(tempconv[0]);
          tempconv = inStructArrayGlo[0]["Lm_ref"];
          scalar_pars.push_back(tempconv[0]);

          // feebacks
          matlab::data::TypedArray<double> feedb = inStructArrayGlo[0]["feedb"];
          std::vector<double> feedbacks(feedb.begin(), feedb.end());
          vector_pars.push_back(feedbacks);

          // modes of action
          matlab::data::TypedArray<double> moac = inStructArrayGlo[0]["moa"];
          std::vector<double> moa(moac.begin(), moac.end());
          vector_pars.push_back(moa);

          // get here some global settings to avoid too much work in read_scen
          matlab::data::TypedArray<double> glo_int_scen = inStructArrayGlo[0]["int_scen"]; // should have just doubles inside 
          std::vector<double> glo_int_scen_vec(glo_int_scen.begin(), glo_int_scen.end());

          matlab::data::TypedArray<double> glo_int_type = inStructArrayGlo[0]["int_type"];
          std::vector<double> glo_int_type_vec(glo_int_type.begin(), glo_int_type.end());


          matlab::data::TypedArray<double> glo_timevar = inStructArrayGlo[0]["timevar"];   // this is also just an array of doubles ([v1, v2])
          std::vector<double> timevar(glo_timevar.begin(), glo_timevar.end());

          matlab::data::TypedArray<double> glo_mf = inStructArrayGlo[0]["MF"];             
          double MF = glo_mf[0];
          scalar_pars.push_back(MF);
		  
		  matlab::data::TypedArray<double> par_a = inStructArrayPar[0]["a"];  // % Weibull background hazard coefficient (-)
		  double a = par_a[0];
		  scalar_pars.push_back(a);
		  
		  // define here the vector int_coll that will be passed afterwards to the solver
          matlab::data::TypedArray<matlab::data::Array> glo_int_coll = inStructArrayGlo[0]["int_coll"];

          auto it = find(glo_int_scen_vec.begin(),glo_int_scen_vec.end(), conc);
		  int int_loc = it - glo_int_scen_vec.begin();
          
          int size_t_int_coll = 0;
          int size_t_int_coll2 = 0;
          //displayOnMATLAB(stream);
          if (timevar[0] == 0){
			  // trick to avoid index error on int_coll calls..
			  // should not affect things
              int_loc = 0;    
          }
          matlab::data::TypedArray<double> int_coll=glo_int_coll[int_loc];
          size_t_int_coll = (int)int_coll.getDimensions()[0];
          size_t_int_coll2 = (int)int_coll.getDimensions()[1];
		  vector <double> int_coll_times;
          int ii =0;
		  while (ii < size_t_int_coll){
			  int_coll_times.push_back(int_coll[ii][0]);
              ii++;
		  }

          int int_type = glo_int_type[int_loc];

          // pass the value to the initial conditions
          //[ state_initialization
          state_type x(4); // in DEB there are 4 states		  
		  // initial conditions
          x[0] = init_states[0]; 
          x[1] = init_states[1];
		  x[2] = init_states[2];
		  x[3] = init_states[3];
          //]

          vector<state_type> x_vec; // states
          vector<double> times;     // times
          
          // Define the stepper type (in this case a dense stepper)
          typedef runge_kutta_dopri5<state_type> stepper_type;

          // CHANGE HERE THE TOLERANCES according to what is in call_deri.m
          double abs_err = AbsErr , rel_err = RelErr , a_x = 1.0 , a_dxdt = 1.0;
		  double max_step = MaxStep;

          // solve the ODE using the stepper already defined. The times are those passed
          // by the user
          size_t steps = integrate_times(make_dense_output(abs_err , rel_err, max_step, stepper_type() ),
                                         DEBderi(scalar_pars,
                                                 vector_pars,
                                                 conc,
                                                 int_coll,
												 int_coll_times,
                                                 int_type,
                                                 timevar,
                                                 matlabPtr2),
                                         x, time_vector, dt, 
                                         push_back_state_and_time( x_vec , times ));
          
          // initialize the arrays to store the output
          matlab::data::TypedArray<double> doubleArray = factory.createArray(
              {times.size(),1}, times.data(), times.data()+times.size());
          matlab::data::TypedArray<double> doubleArray2 = factory.createArray({x_vec.size(),x_vec[0].size()},
                                                                               x_vec[0].data(),
                                                                               x_vec[0].data()+x_vec[0].size());
          
          /* output */
          // manually copy the values of the 4 states (can be improved with another cycle)
          for( int i=0; i<=x_vec.size()-1; i++ )
          {
              // can I avoid this and copy directly the content of the vector?
              doubleArray2[i][0]=x_vec[i][0];
              doubleArray2[i][1]=x_vec[i][1];
              doubleArray2[i][2]=x_vec[i][2];
              doubleArray2[i][3]=x_vec[i][3];
          }
          outputs[0] = doubleArray;  // vector of times
          outputs[1] = doubleArray2; // vector of states
       }
};