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
#ifndef CTBNRLE_FORWARDSAMPLER_H
#define CTBNRLE_FORWARDSAMPLER_H

#include "fwbase.h"
#include "sampler.h"
#include "markovdyn.h"
#include "ctbn.h"
#include "trajectory.h"
#include "structure.h"
#include <vector> 



namespace ctbn {

class ForwardSampler : public FwBase{
public:
    //The caller of the function owns all the pointers
	// but they must remain valie for the lifetime of this object
	ForwardSampler(const Process *pr, const Trajectory *traj, const VarSample *m);
	virtual ~ForwardSampler();

	void SampleTrajectories(std::vector<Trajectory> &tr, std::vector<double> &w,
			int numsamples, Random &rand=randomizer) ;
	ForwardSampler *Clone() const;
	void DiscardBnWeight() {discard_bn_weight = true;}

protected:
	double SampleSingleTrajectory(Trajectory &traj, bool logscale=1, Random &rand=randomizer);
	double SampleDyn(Trajectory &traj, Random &rand=randomizer);
	void CleanCache();
	void MakeCache(const Instantiation &inst); 

	bool discard_bn_weight;

	std::vector<int> nodeindex; 
	std::vector<int> nodeval;
	//std::vector<int> incindex;
	//std::vector<int> currevidval;
	std::vector<int> nextevidval;
};

} // end of ctbn namespace

#endif
