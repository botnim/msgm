////////////////////////////////////////////////////////////////////////////////////////////////////
// file:	stdafx.h
//
// summary:	Windows precompiled header file. We use it to include some of the heavier OpenGM
//	include files inorder to speed up compilation time

#pragma once

// windows headers. 
#ifdef _MSC_VER
#include "targetver.h"	//windows only (used to define maximum windows version)
#include <tchar.h>
#endif

// runtime headers. 
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <random>
#include <tuple>
#include <memory>
#include <ctime>


// OpenGM headers. 
#include "opengm/opengm.hxx"
#include "opengm/graphicalmodel/graphicalmodel.hxx"
#include "opengm/graphicalmodel/graphicalmodel_factor.hxx"
#include "opengm/graphicalmodel/graphicalmodel_hdf5.hxx"
#include "opengm/functions/pottsg.hxx"
#include "opengm/inference/inference.hxx"
#include "opengm/inference/visitors/visitor.hxx"
#include "opengm/graphicalmodel/space/simplediscretespace.hxx"

// Inference algorithms
#include "opengm/inference/bruteforce.hxx"
#include "opengm/inference/trws/trws_trws.hxx"
#include "opengm/inference/lazyflipper.hxx"
#ifndef _MSC_VER
#include "opengm/inference/astar.hxx"
//#include "opengm/inference/lpcplex_mip.hxx"
//#include "opengm/inference/lpcplex_mip_bool.hxx"
#include "opengm/inference/lpcplex_mip_bool.hxx"
#include "opengm/inference/multicut.hxx"
#include "opengm/inference/alphabetaswap.hxx"
#include "opengm/inference/infandflip.hxx"
#include "opengm/inference/dynamicprogramming.hxx"
#include "opengm/inference/external/fastPD.hxx"
//#include "opengm/inference/external/mrflib.hxx"
#include "opengm/inference/external/trws.hxx"
#include "opengm/inference/external/qpbo.hxx"
#include "opengm/inference/mqpbo.hxx"
#include "opengm/inference/external/gco.hxx"
#include "opengm/inference/dualdecomposition/dualdecomposition_bundle.hxx"
#include "opengm/inference/dualdecomposition/dualdecomposition_subgradient.hxx"
#endif 