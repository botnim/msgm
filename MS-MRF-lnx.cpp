////////////////////////////////////////////////////////////////////////////////////////////////////
// file:	cmdline.cpp
// 
// summary:	Implements the command line tool for testing MS-MRF on small datasets and examples
// insead of OpenGM's.

#include "MSMRF.hxx"
#include "stdafx.h"

using namespace std;
using namespace opengm;
using namespace opengm::external;

//template <class GM> void manualSetup(GM& gm, int seed);
template <class GM> void randomSetup(GM& gm, int nVars, int nLabels, double sigma);
template <class GM> void randomSetup8Conn(GM& gm, int nVars, int nLabels, double sigma);
template <class GM> void randomSetup24Conn(GM& gm, int nGridSz, int nLabels, double sigma);
template <class GM> void randomSetup48Conn(GM& gm, int nGridSz, int nLabels, double sigma);
template <class GM> void randomSetupPotts(GM& gm, int nVars, int nLabels, double sigma, int seed);
template <class GM> void randomSetup48ConnDTF(GM& gm, int nGridSz, int nLabels, double sigma);

////////////////////////////////////////////////////////////////////////////////////////////////////
// <summary> Main entry-point for this application.</summary>
// 
// <remarks> The executable accepts as input an HDF5 file with the problem, a relaxation algorithm
// to use and the output file where to store results. Following is example usage:
// 
// 	MS-MRF xxxx.h5 instance_name ICM output.txt
// 	
// </remarks>
// 
// <returns> Exit code: NULL if nothing went wrong.</returns> 

int main(int argc, char* argv[])
{
	// all models of interest are of valuetype double, operator Adder, target minimizer, one of the
	// following function types, and integer index and label types. 
	typedef meta::TypeListGenerator<
		ExplicitFunction<double, size_t, size_t>,
		PottsFunction<double, size_t, size_t>,
		PottsNFunction<double, size_t, size_t>,
		PottsGFunction<double, size_t, size_t>,
		TruncatedSquaredDifferenceFunction<double, size_t, size_t>,
		TruncatedAbsoluteDifferenceFunction<double, size_t, size_t>
	>::type FunctionTypeList;
	//typedef GraphicalModel<double, Adder, ExplicitFunction<double, size_t, size_t> > GM;
	typedef GraphicalModel<double, Adder, FunctionTypeList > GM;

	// setup sub inference params

	//BPS:
	//typedef TRWS<GM> SubInf;	SubInf::Parameter algParam;
	//algParam.doBPS_ = true;
	//algParam.energyType_ = SubInf::Parameter::EnergyType::TABLES;

	//EXPANSION:
	//typedef GCOLIB<GM> SubInf;	SubInf::Parameter algParam;
	//algParam.inferenceType_ = SubInf::Parameter::InferenceType::EXPANSION;
	//algParam.energyType_ = SubInf::Parameter::EnergyType::VIEW;

	//SWAP
	//typedef GCOLIB<GM> SubInf;	SubInf::Parameter algParam;
	//algParam.inferenceType_ = SubInf::Parameter::InferenceType::SWAP;
	//algParam.energyType_ = SubInf::Parameter::EnergyType::VIEW;

	//BUNDLE_A
	//typedef DDDualVariableBlock2<marray::View<typename GM::ValueType,false>> DualBlockType;
	//typedef DualDecompositionBundle<GM, DynamicProgramming<
	//	typename DualDecompositionBase<GM, DualBlockType>::SubGmType, Minimizer>, 
	//	DualBlockType> SubInf;
	//SubInf::Parameter algParam;
	//algParam.minDualWeight_ = 0.1;
	//algParam.maxDualWeight_ = 10000000000;
	//algParam.noBundle_ = true;
	//algParam.stepsizePrimalDualGapStride_ = true;
	//algParam.stepsizeNormalizedSubgradient_ = true;
	//algParam.useHeuristicStepsize_ = false;
	//algParam.stepsizeStride_ = 0.1;
	//algParam.stepsizeExponent_ = 1;
	//algParam.decompositionId_ = DualDecompositionBaseParameter::SPANNINGTREES;
	//algParam.relativeDualBoundPrecision_ = 0;
	//algParam.minimalRelAccuracy_ = 0.0000000001;

	//BUNDLE_H
	//typedef DDDualVariableBlock2<marray::View<typename GM::ValueType,false>> DualBlockType;
	//typedef DualDecompositionBundle<GM, DynamicProgramming<
	//	typename DualDecompositionBase<GM, DualBlockType>::SubGmType, Minimizer>, 
	//	DualBlockType> SubInf;
	//SubInf::Parameter algParam;
	//algParam.minDualWeight_ = 0.1;
	//algParam.maxDualWeight_ = 10000000000;
	//algParam.noBundle_ = true;
	//algParam.useHeuristicStepsize_ = true;
	//algParam.decompositionId_ = DualDecompositionBaseParameter::SPANNINGTREES;
	//algParam.relativeDualBoundPrecision_ = 0;
	//algParam.minimalRelAccuracy_ = 0.0000000001;

	//SUBGRAD_A
	//typedef DDDualVariableBlock2<marray::View<typename GM::ValueType,false>> DualBlockType;
	//typedef DualDecompositionSubGradient<GM, DynamicProgramming<
	//	typename DualDecompositionBase<GM, DualBlockType>::SubGmType, Minimizer>, 
	//	DualBlockType> SubInf;
	//SubInf::Parameter algParam;
	//algParam.useProjectedAdaptiveStepsize_ = true;
	//algParam.stepsizeStride_= 0.1;
	//algParam.stepsizeExponent_ = 1;
	//algParam.stepsizeExponent_ = 0.01;
	//algParam.decompositionId_= DualDecompositionBaseParameter::SPANNINGTREES;
	//algParam.minimalRelAccuracy_ = 0.0000000001;

	//MCA
	//typedef Multicut<GM,Minimizer> SubInf;	SubInf::Parameter algParam;
	//algParam.numThreads_ = 1;

	//TRWS_LF2
	//typedef InfAndFlip<GM,Minimizer,TRWS<GM>> SubInf;	SubInf::Parameter algParam;
	//algParam.subPara_.energyType_ = TRWS<GM>::Parameter::EnergyType::TABLES;

	//SWAP_QPBO
	typedef AlphaBetaSwap<GM,QPBO<GM>> SubInf;
	SubInf::Parameter algParam;

	//TRWS: TRWS<GM>
	//typedef TRWS<GM> SubInf;	SubInf::Parameter algParam;
	//algParam.doBPS_ = false;
	//algParam.energyType_ = SubInf::Parameter::EnergyType::TABLES;

	//QPBO
	//typedef QPBO<GM> SubInf;
	//SubInf::Parameter algParam;

	//FAST PD
	//typedef FastPD<GM> SubInf;
	//SubInf::Parameter algParam;

	//CPLEX
	//typedef LPCplex<GM,Minimizer> SubInf;	SubInf::Parameter algParam;
	//algParam.numberOfThreads_ = 1;
	//algParam.integerConstraint_ = true;

	//LazyFlipper
	//typedef LazyFlipper<GM> SubInf;
	//SubInf::Parameter algParam;

	//MQPBO
	//typedef MQPBO<GM, Minimizer> SubInf;
	//SubInf::Parameter algParam;


	// setup MSMRF params
	MSMRF<GM, SubInf>::Parameter param;
	//algParam.bExhaustItr_ = false;
	param.subInfParam = algParam;
	param.bDegEnt = true;
	param.bUsePostPred = true;
	param.midScaleAlg = MultiscaleBaseParameter::MID_SCALE_ALG::SAME;
	param.nBins = 20;
	param.nVcycles = 3;
	param.nSubInfMaxIter = 1;
	param.nSubInfFinalMaxIter = 10;
	param.nSolveSize = 2;
	param.bImprove = false;
	param.sImproveAlg = "ILP";
	param.bCmptRelax = false;

	//CPLEX standalone params
	//typedef LPCplex<GM,Minimizer> SubInf;	
	//typename SubInf::Parameter algParam;
	//algParam.numberOfThreads_ = 1;
	//algParam.integerConstraint_ = true;
	//algParam.timeLimit_ = 60;


	// Run the tests

	GM* gm;			// current model
	int nBad = 0;	// how many bad instances so far

	for (int j = 5; j <= 5; j++) {

		char *sSig;
		sSig = new char;
		sprintf(sSig, "%02d", j);

		for (int i = 0; i < 1; i++) {
			gm = new GM();

			//hdf5::load(*gm, "/home/stav/MS-MRF/gmb/models/object-seg/objseg-349.h5", "gm");
			//hdf5::load(*gm, "/home/stav/MS-MRF/gmb/models/dtf-chinesechar/TST_test_0043_72_106.h5", "gm");
			//hdf5::load(*gm, "/home/stav/MS-MRF/gmb/models/inpainting-n4/triplepoint4-plain-ring.h5", "gm");
			//hdf5::load(*gm, "/home/stav/MS-MRF/gmb/models/inpainting-n4/triplepoint4-plain-ring-inverse.h5", "gm");
			//hdf5::load(*gm, "/home/stav/MS-MRF/gmb/models/inpainting-n8/triplepoint4-plain-ring.h5", "gm");
			//hdf5::load(*gm, "/home/stav/MS-MRF/gmb/models/brain/t1_icbm_normal_9mm_pn3_rf20.h5", "gm");
			//hdf5::load(*gm, "/home/stav/MS-MRF/gmb/models/mrf-photomontage/family-gm.h5", "gm");

			//manualSetup<GM>(*gm, i+1000);
			randomSetup<GM>(*gm, 30, 2, j*0.1); //INPUT: model, gridSize, numLabels, sigma

			//SubInf cplex(*gm, algParam);


			/////////////
			//// save model
			//char *sIdx;
			//sIdx = new char;
			//sprintf(sIdx, "%03d", i);
			//string sFileName = "/home/stav/MS-MRF/gmb/models/rand-100grid-2labels-";
			//sFileName += sSig;
			//sFileName += "sig/rand-";
			//sFileName += sIdx;
			//sFileName += ".h5";
			//cout << sFileName;
			//hdf5::save(*gm, sFileName, "gm");
			//delete sIdx;
			////////////

			/////////////
			//delete gm;
			//gm = new GM();
			//hdf5::load(*gm, sFileName, "gm");
			/////////////

			//create an inferer around it
			MSMRF<GM, SubInf>* msMRF = 
				new MSMRF<GM, SubInf>(*gm, param);

			msMRF->infer();

			/*
			vector<typename GM::LabelType> start(100*100);
			for (int j = 0; j < 10; j++) {
			*/

			// do the actionZ
			//clock_t t = clock();
			//msMRF->infer();
			//cplex.infer();
			//t = clock() - t;

			//double dMSMRFEnrgy = msMRF->value();
			//double dCPLEXEnrgy = cplex.value();
			//double dCPLEXbound = cplex.bound();

			/*
			cplex.arg(start);
			cplex.setStartingPoint(start.begin());

			}
			*/

			//delete msMRF;
			delete gm;
		}
		delete sSig;
		cout << "Bad: " << nBad;
	}
	return 0;
}



//#define RAND (floor(norm(generator) * 10) / 10)
#define RAND (floor(distribution(generator) * 10) / 10)
#define CUR (r*nGridSz+c)
#define LEFT_OF_CUR (CUR - 1)
#define TOP_OF_CUR (CUR - nGridSz)

// 8 connected
#define TOP_R_OF_CUR (CUR - nGridSz + 1)
#define TOP_L_OF_CUR (CUR - nGridSz - 1)

// 24 connected
#define L_L_OF_CUR (CUR-2)
#define L_L_U_OF_CUR (CUR - 2 - nGridSz)
#define L_L_U_U_OF_CUR (CUR - 2 - 2*nGridSz)
#define L_U_U_OF_CUR (CUR - 1 - 2*nGridSz) 
#define U_U_OF_CUR (CUR-2*nGridSz)
#define R_U_U_OF_CUR (CUR + 1 - 2*nGridSz)
#define R_R_U_U_OF_CUR (CUR + 2 - 2*nGridSz)
#define R_R_U_OF_CUR (CUR + 2 - nGridSz)

// 48 connected
#define L_L_L_OF_CUR (CUR-3)
#define L_L_L_U_OF_CUR (CUR - 3 - nGridSz)
#define L_L_L_U_U_OF_CUR (CUR - 3 - 2*nGridSz)
#define L_L_L_U_U_U_OF_CUR (CUR - 3 - 3*nGridSz) 
#define L_L_U_U_U_OF_CUR (CUR-2-3*nGridSz)
#define L_U_U_U_OF_CUR (CUR - 1 - 3*nGridSz)
#define U_U_U_OF_CUR (CUR + - 3*nGridSz)
#define R_U_U_U_OF_CUR (CUR + 1 - 3*nGridSz)
#define R_R_U_U_U_OF_CUR (CUR + 2 - 3*nGridSz)
#define R_R_R_U_U_U_OF_CUR (CUR + 3 - 3*nGridSz)
#define R_R_R_U_U_OF_CUR (CUR +3 - 2*nGridSz)
#define R_R_R_U_OF_CUR (CUR + 3 - nGridSz)

#define ADD_UNARY(v_, vdRand) \
	v = v_; \
	gm.addVariable(nLabels); \
	ExplicitFunction<double, size_t, size_t> fUnary##v_(shp, shp + 1); \
for (int i = 0; i < nLabels; i++) \
	fUnary##v_(i) = vdRand[i]; \
	gm.addFactor(gm.addFunction(fUnary##v_), &v, &v + 1)

#define ADD_BINARY(vLeft_, vRight_) \
for (int i = 0; i < nLabels*nLabels; i++) \
	vdRand[i] = sigma * RAND; \
	pv[0] = vLeft_; \
	pv[1] = vRight_; \
	ExplicitFunction<double, size_t, size_t> fBinary##vLeft_##vRight_(shp, shp + 2); \
for (int i = 0; i < nLabels; i++) \
for (int j = 0; j < nLabels; j++) \
	fBinary##vLeft_##vRight_(i, j) = vdRand[nLabels * i + j]; \
	gm.addFactor(gm.addFunction(fBinary##vLeft_##vRight_), pv, pv + 2)

#define ADD_BINARY_POTTS(vLeft_, vRight_, dRand) \
	pv[0] = vLeft_; \
	pv[1] = vRight_; \
	ExplicitFunction<double, size_t, size_t> fBinary##vLeft_##vRight_(shp, shp + 2); \
for (int i = 0; i < nLabels; i++) \
for (int j = 0; j < nLabels; j++) \
if (i == j) \
	fBinary##vLeft_##vRight_(i, j) = 0; \
		else \
		fBinary##vLeft_##vRight_(i, j) = dRand; \
		gm.addFactor(gm.addFunction(fBinary##vLeft_##vRight_), pv, pv + 2)

#define ADD_BINARY_DENSE(vLeft_, vRight_) \
if ((vLeft_ >= 0) && (vLeft_ < nVars) && (vRight_ >= 0) && (vRight_ < nVars)) { \
	for (int i = 0; i < nLabels*nLabels; i++) \
		vdRand[i] = sigma * RAND; \
	pv[0] = vLeft_; \
	pv[1] = vRight_; \
	ExplicitFunction<double, size_t, size_t> fBinary##vLeft_##vRight_(shp, shp + 2); \
	for (int i = 0; i < nLabels; i++) \
		for (int j = 0; j < nLabels; j++) \
			fBinary##vLeft_##vRight_(i, j) = vdRand[nLabels * i + j]; \
	gm.addFactor(gm.addFunction(fBinary##vLeft_##vRight_), pv, pv + 2); \
}


template <class GM> void randomSetup48ConnDTF(GM& gm, int nGridSz, int nLabels, double sigma)
{

	int nVars = nGridSz*nGridSz;
	typename GM::IndexType shp[2] = { nLabels, nLabels };
	typename GM::IndexType v, pv[2];
	random_device rdSeed;
	default_random_engine generator(rdSeed());
	//normal_distribution<double> norm(0, 1);
	uniform_real_distribution<double> distribution(0.0, 1.0);


	srand(time(NULL));
	double dPerDT = 90;

	// strong unary term for per% of nodes
	// and no data term for remaining nodes
	for (int v = 0; v < nVars; v++) {
		vector<double> vdRand(nLabels);
		
		//// assign dataterm or skip via random decision
		//double dRand = rand() % 100;
		//if (dRand < dPerDT){
		//	for (int i = 0; i < nLabels; i++)
		//		vdRand[i] = sigma * RAND;
		//}
		//else {
		//	for (int i = 0; i < nLabels; i++)
		//		vdRand[i] = 0;
		//}
		//ADD_UNARY(v, vdRand);

		int nRow = floor(v / nGridSz);
		int nCol = v % nGridSz;
		if ((nRow >= 40) && (nRow < 60) && (nCol >= 40) && (nCol < 60)) {
				for (int i = 0; i < nLabels; i++)
					vdRand[i] = 0;
		}
		else {
			for (int i = 0; i < nLabels; i++)
				vdRand[i] = sigma * RAND;
		}

		ADD_UNARY(v, vdRand);
	}

	// weaker connectivity on pairwise
	sigma = 1;

	// 24-connected
	vector<double> vdRand(nLabels*nLabels);

	//// 4 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(LEFT_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(TOP_OF_CUR, CUR);
		}
	}

	//// 8 connected
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(TOP_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 1; c++) {
			ADD_BINARY(TOP_R_OF_CUR, CUR);
		}
	}

	// 24 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 1; c++) {
			ADD_BINARY(R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 2; c++) {
			ADD_BINARY(R_R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 2; c++) {
			ADD_BINARY(R_R_U_OF_CUR, CUR);
		}
	}


	// 48 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(L_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 1; c++) {
			ADD_BINARY(R_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 2; c++) {
			ADD_BINARY(R_R_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 3; c++) {
			ADD_BINARY(R_R_R_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 3; c++) {
			ADD_BINARY(R_R_R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 3; c++) {
			ADD_BINARY(R_R_R_U_OF_CUR, CUR);
		}
	}
}


template <class GM> void randomSetup24Conn(GM& gm, int nGridSz, int nLabels, double sigma)
{

	int nVars = nGridSz*nGridSz;
	typename GM::IndexType shp[2] = { nLabels, nLabels };
	typename GM::IndexType v, pv[2];
	random_device rdSeed;
	default_random_engine generator(rdSeed());
	//normal_distribution<double> norm(0, 1);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	for (int v = 0; v < nVars; v++) {
		vector<double> vdRand(nLabels);
		for (int i = 0; i < nLabels; i++)
			vdRand[i] = (1-sigma) * RAND;
		ADD_UNARY(v, vdRand);
	}

	// 24-connected
	vector<double> vdRand(nLabels*nLabels);

	//// 4 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(LEFT_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(TOP_OF_CUR, CUR);
		}
	}

	//// 8 connected
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(TOP_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 1; c++) {
			ADD_BINARY(TOP_R_OF_CUR, CUR);
		}
	}

	// 24 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 1; c++) {
			ADD_BINARY(R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 2; c++) {
			ADD_BINARY(R_R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 2; c++) {
			ADD_BINARY(R_R_U_OF_CUR, CUR);
		}
	}
}


template <class GM> void randomSetup48Conn(GM& gm, int nGridSz, int nLabels, double sigma)
{

	int nVars = nGridSz*nGridSz;
	typename GM::IndexType shp[2] = { nLabels, nLabels };
	typename GM::IndexType v, pv[2];
	random_device rdSeed;
	default_random_engine generator(rdSeed());
	//normal_distribution<double> norm(0, 1);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	for (int v = 0; v < nVars; v++) {
		vector<double> vdRand(nLabels);
		for (int i = 0; i < nLabels; i++)
			vdRand[i] = (1-sigma) * RAND;
		ADD_UNARY(v, vdRand);
	}

	// 48-connected
	vector<double> vdRand(nLabels*nLabels);

	//// 4 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(LEFT_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(TOP_OF_CUR, CUR);
		}
	}

	//// 8 connected
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(TOP_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 1; c++) {
			ADD_BINARY(TOP_R_OF_CUR, CUR);
		}
	}

	// 24 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 1; c++) {
			ADD_BINARY(R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 2; c++) {
			ADD_BINARY(R_R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 2; c++) {
			ADD_BINARY(R_R_U_OF_CUR, CUR);
		}
	}

	// 48 connected
	for (int r = 0; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 3; c < nGridSz; c++) {
			ADD_BINARY(L_L_L_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 2; c < nGridSz; c++) {
			ADD_BINARY(L_L_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 1; c < nGridSz; c++) {
			ADD_BINARY(L_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz; c++) {
			ADD_BINARY(U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz-1; c++) {
			ADD_BINARY(R_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz-2; c++) {
			ADD_BINARY(R_R_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 3; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 3; c++) {
			ADD_BINARY(R_R_R_U_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 2; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 3; c++) {
			ADD_BINARY(R_R_R_U_U_OF_CUR, CUR);
		}
	}
	for (int r = 1; r < nGridSz; r++) {
		for (int c = 0; c < nGridSz - 3; c++) {
			ADD_BINARY(R_R_R_U_OF_CUR, CUR);
		}
	}
}


template <class GM> void randomSetup(GM& gm, int nGridSz, int nLabels, double sigma)
{
	int nVars = nGridSz*nGridSz;
	typename GM::IndexType shp[2] = { nLabels, nLabels };
	typename GM::IndexType v, pv[2];
	random_device rdSeed;
	default_random_engine generator(rdSeed());
	//generator.seed(rdSeed());
	//normal_distribution<double> norm(0, 1);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	for (int v = 0; v < nVars; v++) {
		vector<double> vdRand(nLabels);
		for (int i = 0; i < nLabels; i++)
			vdRand[i] = (1-sigma) * RAND;
		ADD_UNARY(v, vdRand);
	}

	for (int r = 0; r < nGridSz; r++) {
		if (r == 0) {
			for (int c = 1; c < nGridSz; c++) {
				vector<double> vdRand(nLabels*nLabels);
				//for (int i = 0; i < nLabels*nLabels; i++)
				//	vdRand[i] = sigma * RAND;
				ADD_BINARY(LEFT_OF_CUR, CUR);
			}
		}
		else {
			for (int c = 0; c < nGridSz; c++) {
				vector<double> vdRand(nLabels*nLabels);
				//for (int i = 0; i < nLabels*nLabels; i++)
					//vdRand[i] = sigma * RAND;
				if (c > 0) {
					ADD_BINARY(LEFT_OF_CUR, CUR);
				}

				ADD_BINARY(TOP_OF_CUR, CUR);
			}
		}
	}
}

template <class GM> void randomSetupPotts(GM& gm, int nGridSz, int nLabels, double sigma, int seed)
{
	int nVars = nGridSz*nGridSz;
	typename GM::IndexType shp[2] = { nLabels, nLabels };
	typename GM::IndexType v, pv[2];
	default_random_engine generator(1);
	generator.seed(seed);
	//normal_distribution<double> norm(0, 1);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	for (int v = 0; v < nVars; v++) {
		vector<double> vdRand(nLabels);
		for (int i = 0; i < nLabels; i++)
			vdRand[i] = RAND;
		ADD_UNARY(v, vdRand);
	}

	for (int r = 0; r < nGridSz; r++) {
		if (r == 0) {
			for (int c = 1; c < nGridSz; c++) {
				double dRand;
				dRand = sigma * RAND;
				ADD_BINARY_POTTS(LEFT_OF_CUR, CUR, dRand);
			}
		}
		else {
			for (int c = 0; c < nGridSz; c++) {
				double dRand;
				dRand = sigma * RAND;
				if (c > 0) {
					ADD_BINARY_POTTS(LEFT_OF_CUR, CUR, dRand);
				}

				ADD_BINARY_POTTS(TOP_OF_CUR, CUR, dRand);
			}
		}
	}
}


template <class GM> void randomSetup8Conn(GM& gm, int nGridSz, int nLabels, double sigma)
{
	int nVars = nGridSz*nGridSz;
	typename GM::IndexType shp[2] = { nLabels, nLabels };
	typename GM::IndexType v, pv[2];
	random_device rdSeed;
	default_random_engine generator(rdSeed());
	//normal_distribution<double> norm(0, 1);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	for (int v = 0; v < nVars; v++) {
		vector<double> vdRand(nLabels);
		for (int i = 0; i < nLabels; i++)
			vdRand[i] = (1-sigma) * RAND;
		ADD_UNARY(v, vdRand);
	}

	vector<double> vdRand(nLabels*nLabels);
	for (int r = 0; r < nGridSz; r++) {
		if (r == 0) {
			for (int c = 1; c < nGridSz; c++) {
				ADD_BINARY(LEFT_OF_CUR, CUR);
			}
		}
		else {
			for (int c = 0; c < nGridSz; c++) {
				
				// connect left variable
				if (c > 0) {
					ADD_BINARY(LEFT_OF_CUR, CUR);
				}

				// connect top variable
				ADD_BINARY(TOP_OF_CUR, CUR);

				// connect top right variable
				if (c < nGridSz - 1) {
					ADD_BINARY(TOP_R_OF_CUR, CUR);
				}

				// connect top left variable
				if (c > 0) {
					ADD_BINARY(TOP_L_OF_CUR, CUR);
				}

			}
		}
	}
}
