/* 
 * File:   Source file for class Solution
 * Author: A. Herrán
 */

#include "Solution.h"

//------------------------
// Constructor/destructor
//------------------------
Solution::Solution()  { };
Solution::~Solution() { };

//----------------------------------------------------------------------
// Construct Operators: Random, GreedyRandom, RandomGreedy
//----------------------------------------------------------------------
void Solution::Construct(double alpha) {

 // Initializes parameters structures
    vecUns cluster;
    sol = matUns(instance.nK, cluster);
    obj = 0.0;

 // Constructs a new solution
    if     (parameters.CONm=="random"  ) ConstructiveRND();
    else if(parameters.CONm=="greedy-a") ConstructiveGRD(alpha);
    else if(parameters.CONm=="greedy-r") ConstructiveGRD(alpha);
    else {
        cout << "\n Wrong construction method.";
        cin.get();
        exit(1);
    }                 
    
}

//-------------------------------------------------------------------------------
// Initializes parameters structures needed for an efficient evaluation of moves
//-------------------------------------------------------------------------------
void Solution::DataInit() {

    delta = matDbl(instance.nV,vecDbl(instance.nK, 0.0));
    for(unsigned i=0; i<instance.nV; i++) {
        for(unsigned k=0; k<instance.nK; k++) {
            for(unsigned j: sol[k]) delta[i][k] += instance.D[i][j] + instance.D[j][i];
        }
    }

    sigma = vecDbl(instance.nK, 0.0);
    w     = matDbl(instance.nK, vecDbl(instance.nT, 0.0));
    for(unsigned k=0; k<instance.nK; k++) {
        for(unsigned t=0; t<instance.nT; t++) {
            for(unsigned i: sol[k]) w[k][t] += instance.W[i][t];
            if(w[k][t] < instance.WL[k][t]) sigma[k] += instance.WL[k][t] - w[k][t];
            if(w[k][t] > instance.WU[k][t]) sigma[k] += w[k][t] - instance.WU[k][t];
        }
    }

}

//----------------------------------------------------------------------
// Eval Operators: Full solution and full generic move
//----------------------------------------------------------------------

void Solution::Evaluate(bool objective) {

    if(parameters.LSe) DataInit();

    obj = 0.0;
    if(objective) {
        for(unsigned k=0; k<instance.nK; k++) {
            for(unsigned i=0; i<sol[k].size()-1; i++) {
                for(unsigned j=i+1; j<sol[k].size(); j++) {
                    obj += instance.D[sol[k][i]][sol[k][j]] + instance.D[sol[k][j]][sol[k][i]];
                }
            }
        }
    }
    else {
        if(parameters.LSe) for(unsigned k=0; k<instance.nK; k++) obj += sigma[k];
        else {
            for(unsigned k=0; k<instance.nK; k++) {
                for(unsigned t=0; t<instance.nT; t++) {
                    double wkt = 0.0;
                    for(unsigned i: sol[k]) wkt += instance.W[i][t];
                    if(wkt < instance.WL[k][t]) obj += instance.WL[k][t] - wkt;
                    if(wkt > instance.WU[k][t]) obj += wkt - instance.WU[k][t];
                }
            }
        }
    }

}

double Solution::EvaluateMove(bool objective) {

    double objaux = 0.0;
    if(objective) {
        for(unsigned k=0; k<instance.nK; k++) {
            for(unsigned i=0; i<sol[k].size()-1; i++) {
                for(unsigned j=i+1; j<sol[k].size(); j++) {
                    objaux += instance.D[sol[k][i]][sol[k][j]] + instance.D[sol[k][j]][sol[k][i]];
                }
            }
        }
    }
    else {
        for(unsigned k=0; k<instance.nK; k++) {
            for(unsigned t=0; t<instance.nT; t++) {
                double wkt = 0.0;
                for(unsigned i: sol[k]) wkt += instance.W[i][t];
                if(wkt < instance.WL[k][t]) objaux += instance.WL[k][t] - wkt;
                if(wkt > instance.WU[k][t]) objaux += wkt - instance.WU[k][t];
            }
        }
    }

    return objaux - obj;
}

double Solution::EvaluateHMP() {

    double objaux = 0.0;
    for(unsigned k=0; k<instance.nK-1; k++) {
        for(unsigned kp=k+1; kp<instance.nK; kp++) {
            for(unsigned i=0; i<sol[k].size(); i++) {
                for (unsigned ip=0; ip<sol[kp].size(); ip++) {
                    objaux -= instance.D[sol[k][i]][sol[kp][ip]];
                    objaux -= instance.D[sol[kp][ip]][sol[k][i]];
                }
            }
        }
    }

    return objaux;
}

//------------------------------------------------------------------
// Evaluates the feasibility of a move between two clusters
//------------------------------------------------------------------
bool Solution::FeasibleMove(unsigned k1, unsigned k2) {

    bool feasible = true;

    for(unsigned t=0; t<instance.nT; t++) {
        double wkt = 0.0;
        for(unsigned i: sol[k1]) wkt += instance.W[i][t];
        if((wkt < instance.WL[k1][t]) || (wkt > instance.WU[k1][t])) { feasible = false; break; }
        wkt = 0.0;
        for(unsigned i: sol[k2]) wkt += instance.W[i][t];
        if((wkt < instance.WL[k2][t]) || (wkt > instance.WU[k2][t])) { feasible = false; break; }
    }

    return feasible;
}

//----------------------------------------------------------------------
// Neighborhood exploration strategies
//----------------------------------------------------------------------

bool Solution::Explore(bool objective) {

    bool improved;

    if     (parameters.LSm =="ins") improved = ExploreInsert  (objective);
    else if(parameters.LSm =="exc") improved = ExploreExchange(objective);
    else if(parameters.LSm =="ext") improved = ExploreExtended(objective);
    else {
        cout << "\n Wrong neighborhood exploration method.";
        cin.get();
        exit(1);
    }

    return improved;
}
