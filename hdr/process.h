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
#ifndef CTBNRLE_PROCESS_H
#define CTBNRLE_PROCESS_H

#include "dynamics.h"
#include "rv.h"



namespace ctbn {

// Process is a ABC for any general dynamic process
class Process : public StreamObj {
	//SOBJCLASSDECL(Process)
public:
	virtual ~Process();

	// virtual copy constructor
	virtual Process *Clone() const = 0;

	// the amalgamtion operation
	virtual void Mult(const Process *p) = 0;
	
	// returns a suitable sufficient statistics object
	virtual SS *BlankSS() const = 0;
	
	// sets parameters from ML estimate
	virtual void Maximize(const SS *ss) = 0;

    // returns log-likelihood from sufficient statistics
    virtual double LLH(const SS *ss) const = 0;

	// sample a trajectory of events from the process
	virtual void Sample(Trajectory &tr, 
				Random &rand=randomizer) const = 0;

	SERIAL_START_V(Process)
	SERIAL_END
};

} // end of ctbn namespace

#endif

