/* 
 * File:   Source file for class GRASP
 * Author: A. Herrán
 */

#include "../metaheuristics/GRASP.h"

//------------------------
// Constructor/destructor
//------------------------
GRASP::GRASP(): Metaheuristic() {

    solutione.obj = std::numeric_limits<double>::max();

    alpha = parameters.GRASPa;
    m     = parameters.GRASPm;
    block = parameters.GRASPb;
    delta = parameters.GRASPd;

 // Reactive GRASP
    if(parameters.CONm=="greedy-r") {
        m++;
        valpha.emplace_back(0.0);
        for(int i=1; i<m; i++) valpha.emplace_back(valpha[i-1] + 1.0/(double)(m-1));
        n   = vecUns(m,0);
        p   = vecDbl(m,1.0/(double)m);
        sum = vecDbl(m,0.0);
        q   = vecDbl(m,0.0);
    }

}

GRASP::~GRASP() = default;

//----------------------------------------------
// Performs one iteration of the meta-heuristic
//----------------------------------------------
void GRASP::Iterate() {

 // Reactive GRASP
    int i;
    if(parameters.CONm=="greedy-r") {
        if (k < m) {
            i = k;
        } else {
            vecDbl roulette = vecDbl(m, 0.0);
            roulette[0] = p[0];
            for(int j=1; j<m; j++) roulette[j] = roulette[j-1] + p[j];
            double shot = rndDbl(0.0, 1.0);
            i = 0;
            while (roulette[i] < shot) i++;
        }
        n[i]++;
        alpha = valpha[i];
    }

 // First phase of the Local Search procedure
    bool feasible = false;
    while(!feasible) {
        feasible = true;
        solutiona.Construct(alpha);
        solutiona.Evaluate(false);
        while(solutiona.Explore(false)) continue;
        if(myRound(solutiona.obj)>0.0) feasible = false;
    }   
    solutiona.Evaluate(true);

 // Second phase of the Local Search procedure
    while(solutiona.Explore(true)) continue;

 // Update incumbent solution
    if(solutiona<solutione) solutione = solutiona;

 // Reactive GRASP
    if(parameters.CONm=="greedy-r") {
        sum[i] += solutiona.obj;
        if((k+1)%block==0) {
            double sumq = 0.0;
            for(int i=0; i<m; i++) {
                double  Ai = sum[i]/(double)n[i];
                q[i] = pow(solutione.obj/Ai,delta);
                sumq += q[i];
            }
            for(int i=0; i<m; i++) p[i] = q[i]/sumq;
        }
    }

    k++;

}
