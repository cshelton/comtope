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
#ifndef CTBNRLE_NONEXPSUFFSTATSQUERY_H
#define CTBNRLE_NONEXPSUFFSTATSQUERY_H


// A class for performing queries for sufficient statistics where we have 
// complete trajectories.

#include "suffstatsquery.h"
#include "trajectory.h"
#include <vector>

namespace ctbn {

class NonExpSuffStatsQuery : public SuffStatsQuery {
public:
	NonExpSuffStatsQuery();
	virtual ~NonExpSuffStatsQuery();
	virtual SS* GetP0SS();
	virtual SS* GetDSS();
};

} // end of ctbn namespace

#endif
