////////////////////////////////////////////////////////////////////////////////////////////////////
// file:	cmdline.cpp
// 
// summary:	Implements the command line tool for testing MS-MRF on small datasets and examples
// insead of OpenGM's.

#include "stdafx.h"
#include "MSMRF.hxx"

using namespace std;
using namespace opengm;

template <class GM> void manualSetup(GM& gm, int seed);
template <class GM> void randomSetup(GM& gm, int nVars, int seed);

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
	/*if (argc != 4) {
		cerr << "MS-MRF command line takes 4 arguments.";
		return -1;
	}*/

	
	// all models of interest are of valuetype double, operator Adder, target minimizer, one of the
	// following function types, and integer index and label types. 
	/*typedef meta::TypeListGenerator<
      ExplicitFunction<double, size_t, size_t>,
      PottsFunction<double, size_t, size_t>,
      PottsNFunction<double, size_t, size_t>,
      PottsGFunction<double, size_t, size_t>,
      TruncatedSquaredDifferenceFunction<double, size_t, size_t>,
      TruncatedAbsoluteDifferenceFunction<double, size_t, size_t> 
      >::type FunctionTypeList;*/
	
	typedef GraphicalModel<double, Adder, ExplicitFunction<double, size_t, size_t> > Model;

	//QPBO
	typedef Bruteforce<Model, Minimizer> SubInf;
	SubInf::Parameter algParam;

	// setup MSMRF params
	MSMRF<Model, SubInf>::Parameter param;
	param.bDegEnt = false;
	param.bUsePostPred = false;
	param.subInfParam = algParam;
	param.midScaleAlg = MultiscaleBaseParameter::MID_SCALE_ALG::NONE;
	param.nBins = 20;
	param.nVcycles = 1;
	param.nSubInfMaxIter = 1;
	param.nSubInfFinalMaxIter = 100;
	param.nSolveSize = 100;
	
	// load the data model. 
	Model* gm;

	int nBad = 0;
	for (int i = 0; i < 2000; i++) {
		gm = new Model();

		//hdf5::load(*gm, "S:\\Dropbox\\Projects\\ms-mrf\\ms-mrf\\datasets\\chinese\\TST_test_0009_64_80.h5", "gm");
		//size_t n = gm.numberOfFactors();
		//manualSetup<Model>(*gm, i+2);
		randomSetup<Model>(*gm, 30, 1);

		////create an inferer around it
		//MSMRF<Model, SubInf>* msMRF = 
		//	new MSMRF<Model, SubInf>(*gm, param);

		//// do the action
		//size_t t = GetTickCount();
		//msMRF->infer();
		//t = GetTickCount() - t;
		////cout << t;
		//double dMSMRFEnrgy = msMRF->value();

		//// solve using brute
		//Bruteforce<Model, Minimizer>* brute = new Bruteforce<Model, Minimizer>(*gm);
		//brute->infer();
		//double dBruteEnrgy = brute->value();
	
		//if (floor(dBruteEnrgy*1e6 + 0.5) != floor(dMSMRFEnrgy*1e6 + 0.5)) {
		//	assert(floor(dBruteEnrgy*1e6 + 0.5) < floor(dMSMRFEnrgy*1e6 + 0.5));
		//	//cout << "----------" << endl;
		//	//cout << dMSMRFEnrgy << endl;
		//	//cout << dBruteEnrgy << endl;
		//	cout << "index: " << i << ", diff: " << dMSMRFEnrgy - dBruteEnrgy << endl;
		//	nBad++;
		//}

		delete msMRF;
		delete gm;
	}

	// output the argmin
	/*vector<size_t> argmin;
	msMRF.arg(argmin);
	for(size_t variable = 0; variable < gm.numberOfVariables(); ++variable) {
		cout << "x" << variable << "=" << argmin[variable] << "\n";
	} 
*/
	cout << "Bad: " << nBad;
	return 0;
}


#define ADD_UNARY(v_, val0, val1) \
	v = v_; \
	gm.addVariable(nLabels); \
	ExplicitFunction<double, size_t, size_t> fUnary##v_(shp, shp + 1); \
	fUnary##v_(0) = val0; \
	fUnary##v_(1) = val1; \
	gm.addFactor(gm.addFunction(fUnary##v_), &v, &v + 1)

#define ADD_BINARY(vLeft_, vRight_, val00, val01, val10, val11) \
	pv[0] = vLeft_; \
	pv[1] = vRight_; \
	ExplicitFunction<double, size_t, size_t> fBinary##vLeft_##vRight_(shp, shp + 2); \
	fBinary##vLeft_##vRight_(0, 0) = val00; \
	fBinary##vLeft_##vRight_(0, 1) = val01; \
	fBinary##vLeft_##vRight_(1, 0) = val10; \
	fBinary##vLeft_##vRight_(1, 1) = val11; \
	gm.addFactor(gm.addFunction(fBinary##vLeft_##vRight_), pv, pv + 2)

#define ADD_UNARY_3(v_, val0, val1, val2) \
	v = v_; \
	gm.addVariable(nLabels); \
	ExplicitFunction<double, size_t, size_t> fUnary##v_(shp, shp + 1); \
	fUnary##v_(0) = val0; \
	fUnary##v_(1) = val1; \
	fUnary##v_(2) = val2; \
	gm.addFactor(gm.addFunction(fUnary##v_), &v, &v + 1)

#define ADD_BINARY_3(vLeft_, vRight_, i_, vdRandVals_)\
	pv[0] = vLeft_; \
	pv[1] = vRight_; \
	ExplicitFunction<double, size_t, size_t> fBinary##vLeft_##vRight_(shp, shp + 2); \
	fBinary##vLeft_##vRight_(0, 0) = vdRandVals_[12 + 9* i_]; \
	fBinary##vLeft_##vRight_(0, 1) = vdRandVals_[12 + 9 * i_ + 1]; \
	fBinary##vLeft_##vRight_(0, 2) = vdRandVals_[12 + 9 * i_ + 2]; \
	fBinary##vLeft_##vRight_(1, 0) = vdRandVals_[12 + 9 * i_ + 3]; \
	fBinary##vLeft_##vRight_(1, 1) = vdRandVals_[12 + 9 * i_ + 4]; \
	fBinary##vLeft_##vRight_(1, 2) = vdRandVals_[12 + 9 * i_ + 5]; \
	fBinary##vLeft_##vRight_(2, 0) = vdRandVals_[12 + 9 * i_ + 6]; \
	fBinary##vLeft_##vRight_(2, 1) = vdRandVals_[12 + 9 * i_ + 7]; \
	fBinary##vLeft_##vRight_(2, 2) = vdRandVals_[12 + 9 * i_ + 8]; \
	gm.addFactor(gm.addFunction(fBinary##vLeft_##vRight_), pv, pv + 2)

template <class GM> void manualSetup(GM& gm, int seed)
{

	GM::IndexType nVars = 4, nLabels = 2;
	GM::IndexType shp[2] = {nLabels, nLabels};
	GM::IndexType v, pv[2];

	/*ADD_UNARY(0, -0.6,  0.6);
	ADD_UNARY(1,  1.2, -0.9);
	ADD_UNARY(2, -0.8,  0.8);
	ADD_UNARY(3, -1.1, -1.5);
	ADD_UNARY(4, -0.8,  0.9);
	ADD_UNARY(5, -0.6, -0.2);
	ADD_UNARY(6, -0.6,  0.2);
	ADD_UNARY(7,  0.2, -2  );
	ADD_UNARY(8, -0.2, -1.3);

	ADD_BINARY(0,1,   1.2, -0.3,    2,  0.6);
	ADD_BINARY(1,2,   1.8, -1.9,  1.8, -1.1);
	ADD_BINARY(3,4,  -0.4, -1.4,  1.4, -0.3);
	ADD_BINARY(4,5,   1.3,  0.5,  0.1, -0.3);
	ADD_BINARY(6,7,   0.8, -0.5,  0.9,  1.8);
	ADD_BINARY(7,8,   0.6,  0.6, -0.6, -0.5);
	ADD_BINARY(0,3,  -0.2, -1  ,  0.6, -0.3);
	ADD_BINARY(1,4,  -1.2,  0.9, -0.7, -1.5);
	ADD_BINARY(2,5,  -0.6, -1.5,    0,  0.3);
	ADD_BINARY(3,6,  -2  ,  0.4,  0.2, -1.3);
	ADD_BINARY(4,7,  -0.5, -0.3,  0.6, -1.9);
	ADD_BINARY(5,8,   0  ,  0.1,  0.3,  0.4);

	return;*/

	//TODO: find how to init
	default_random_engine generator(1);
	generator.seed(seed);
	normal_distribution<double> norm(0,1);
	
	if (nLabels == 2) {
		// 2 labels

		// fill vector with random values
		vector<double> vdRandVals(48);
		for (int i = 0; i < 48; i++) {
			vdRandVals[i] = floor(norm(generator) * 10) / 10;
		}

		ADD_UNARY(0, vdRandVals[0], vdRandVals[1]);
		ADD_UNARY(1, vdRandVals[2], vdRandVals[3]);
		ADD_UNARY(2, vdRandVals[4], vdRandVals[5]);
		ADD_UNARY(3, vdRandVals[6], vdRandVals[7]);

		ADD_BINARY(0, 1, vdRandVals[8], vdRandVals[9], vdRandVals[10], vdRandVals[11]);
		ADD_BINARY(2, 3, vdRandVals[12], vdRandVals[13], vdRandVals[14], vdRandVals[15]);
		ADD_BINARY(0, 2, vdRandVals[16], vdRandVals[17], vdRandVals[18], vdRandVals[19]);
		ADD_BINARY(1, 3, vdRandVals[20], vdRandVals[21], vdRandVals[22], vdRandVals[23]);
	}
	else {
		//// 3 labels

		// fill vector with random values
		vector<double> vdRandVals(48);
		for (int i = 0; i < 48; i++) {
			vdRandVals[i] = floor(norm(generator) * 10) / 10;
		}

		// construct model
		for (int i = 0; i < 4; i++) {
			ADD_UNARY_3(i, vdRandVals[3 * i], vdRandVals[3 * i + 1], vdRandVals[3 * i + 2]);
			//cout << vdRandVals[3 * i] << "\t"
			//	<< vdRandVals[3 * i + 1] << "\t"
			//	<< vdRandVals[3 * i + 2] << endl;
		}

		//cout << "--------" << endl;

		for (int i = 0; i < 4; i++) {
			if (i < 3) {
				int j = i + 1;
				ADD_BINARY_3(i, j, i, vdRandVals);
				//for (int j = 12; j < 21; j++){
				//	cout << vdRandVals[j] << "\t";
				//}
			}
			else {
				ADD_BINARY_3(0, 3, i, vdRandVals);
			}		
		}
	}

}

#define RAND (floor(norm(generator) * 10) / 10)
#define CUR (r*nGridSz+c)
#define LEFT_OF_CUR (CUR - 1)
#define TOP_OF_CUR (CUR - nGridSz)

template <class GM> void randomSetup(GM& gm, int nGridSz, int seed)
{
	int nVars = nGridSz*nGridSz;
	int nLabels = 2;
	GM::IndexType shp[2] = {nLabels, nLabels};
	GM::IndexType v, pv[2];
	default_random_engine generator(1);
	generator.seed(seed);
	normal_distribution<double> norm(0,1);

	for (int v = 0; v < nVars; v++) {
		ADD_UNARY(v, RAND, RAND);
	}

	for (int r = 0; r < nGridSz; r++) {
		if (r == 0) {
			for (int c = 1; c < nGridSz; c++) {
				ADD_BINARY(LEFT_OF_CUR, CUR, RAND, RAND, RAND, RAND);
			}
		}
		else {
			for (int c = 0; c < nGridSz; c++) {
				if (c > 0) {
					ADD_BINARY(LEFT_OF_CUR, CUR, RAND, RAND, RAND, RAND);
				}
				
				ADD_BINARY(TOP_OF_CUR, CUR, RAND, RAND, RAND, RAND);
			}
		}
	}
}