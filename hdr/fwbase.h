/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * Authored by Yu Fan, William Lam, Joon Lee, Christian Shelton, and Jing Xu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CTBNRLE_FWBASE_H
#define CTBNRLE_FWBASE_H

#include "sampler.h"
#include "markovdyn.h"
#include "ctbn.h"
#include "trajectory.h"
#include "structure.h"
#include <vector> 



namespace ctbn {

//FwBase is for forward sampling
class FwBase : public Sampler{
public:
	virtual ~FwBase(){};
	virtual void SetProcess(const Process *pr);
	virtual void SampleTrajectories(std::vector<Trajectory> &tr, std::vector<double> &w, 
            int numsamples, Random &rand=randomizer) = 0;
protected:

	virtual double SampleInitial(Instantiation &initval, Random &rand=randomizer);
	virtual double SampleDyn(Trajectory &tr, Random &rand=randomizer) = 0;
	//sample the next transition time for a variable
	virtual double SampleTime(const Instantiation &curri, int id, 
			int parindex, int nodeval,
			double t, Random &rand=randomizer) ;
	//weight contribution of sampling transition time
	virtual double SampleTimeWeight(const Instantiation &curri, int nodeid, int transitionid, 
			int parindex, int nodeval,
			double currt, double nextt);
	//sample the next state for a variable and 
	//update the corresponding weight contribution
	virtual int SampleTransition(Instantiation &inst, int id,
			int parindex, int nodeval,
			double t, double &weight, Random &rand=randomizer);
	//weight contribution of sampling the next state
	/*virtual double SampleTransitionWeight(const Instantiation &nextinst, int transitionid, int oldval, 
			int parindex, int nodeval,
			int evidencetype, Trajectory::Index &nextindex, int e,
			double currt, double nextt);*/
	virtual void UndefineTime(int id, std::vector<int> &undefinelist);

	Structure str;  
	std::vector<std::vector<int> > markovblanket;
	std::vector<const MarkovDyn*> nodes;
	//std::vector<Trajectory::Index> evidindex;
	bool forceflip;
	std::vector<bool> changeindex;
};

} // end of ctbn namespace

#endif
