////////////////////////////////////////////////////////////////////////////////////////////////////
// file:	MSMRF.hxx
// 
// summary:	Definition and implementation of the MSMRF inference algorithm. 

#pragma once
#include "stdafx.h"
#include "MultSclGM.hxx"

// constants
#define COARSE_VERTEX_RATIO	0.7		// the ratio between the number of vertices in the coarse model to fine
#define TRUE 1
#define FALSE 0

using namespace std;
#ifndef _MSC_VER
using namespace opengm::external;
#endif

namespace opengm {

struct MultiscaleBaseParameter {
	enum MID_SCALE_ALG {
		NONE,
		SAME,
		SWAP,
		ICM
	};
	
	bool bCmptRelax;		///< Flag for interpolation via compatible relaxations
	bool bImprove;			///< Flag to improve solution from 'BPS' or 'LBP-LF'
	string sImproveAlg;		///< improving solution from 'BPS' or 'LBP-LF'
	bool bDegEnt;			///< Flag for degree normalization
	bool bUsePostPred;		///< Flag for using post predictive
	int nBins;				///< The number of bins (0 for no binning)
	int nVcycles;			///< Number of v-cycles
	int nSubInfMaxIter;		///< Maximum iterations for sub-inference algorithm
	int nSubInfFinalMaxIter;	///< Max iterations at the finest level
	MID_SCALE_ALG midScaleAlg;  ///< algorithm to use for mid-scales
	int nSolveSize;

	MultiscaleBaseParameter() :
		bDegEnt(false), nBins(20), bUsePostPred(false), nVcycles(1),
		nSubInfMaxIter(5), nSubInfFinalMaxIter(5), midScaleAlg(SWAP),
		nSolveSize(2), bImprove(false), bCmptRelax(false) {}
};

// <summary> MSMRF Class - MSMRF inference algorithm.</summary>
// <remarks> This class inherits from OpenGM's standard Inference interface and implements the
//  multiscale inference algorithm as suggested in:
//  @@@@@@@@@.</remarks>
// <typeparam name="GM"> 	Type of the graphical model to be processed (only min-sum models are
// 							supported).</typeparam>
// <typeparam name="ALG"> Type of the single scale algorithm.</typeparam>
template<class GM, class ALG>
class MSMRF :
	public Inference<GM, Minimizer>
{
public:
	// get the OpenGM typedefs
	typedef GM GraphicalModelType;
	OPENGM_GM_TYPE_TYPEDEFS;
	typedef VerboseVisitor<MSMRF<GM,ALG>> VerboseVisitorType;
    typedef EmptyVisitor<MSMRF<GM,ALG>> EmptyVisitorType;
    typedef TimingVisitor<MSMRF<GM,ALG>> TimingVisitorType;

	typedef MultSclGM<GM> MSGM;
	typedef typename MSGM::Vertex Vertex;
	typedef typename MSGM::Edge Edge;
	typedef typename MSGM::ExFunction ExFunction;

	typedef typename ALG::Parameter SubParam;

	struct Parameter : public MultiscaleBaseParameter {
		
		SubParam subInfParam;

		Parameter(MultiscaleBaseParameter& param) :
			MultiscaleBaseParameter(param) {}

		Parameter() : MultiscaleBaseParameter() {}
		
	};

	// <summary> Constructor for graphical model.</summary>
	// <remarks> Constructs an MS-MRF inference algorithm for the specified GraphicalModel.</remarks>
	// <param name="gm"> The GraphicalModel to operate on.</param>
	MSMRF(GM& gm, const Parameter& param = Parameter());

	// <summary> Destructor.</summary>
	~MSMRF(void);

	
	// Inference interface implementation
	std::string name() const { return "MS-MRF"; }
	const GM& graphicalModel() const { return gmFine; }

	// <summary> Infers the solution to the graphical model.</summary>
	// <remarks> This follows the guidelines of the Inference interface.</remarks>
	// <returns> An InferenceTermination status.</returns>
	InferenceTermination infer() { EmptyVisitorType v; return infer(v); };
	template<class VisitorType> InferenceTermination infer(VisitorType&);

	// <summary> The infered argmin.</summary>
	// <remarks> Returns the state of the variables yielding the minimal value.</remarks>
	// <param name="vLabels"> [in,out] Vector of labels (indexed according to variable indices).</param>
	// <param name="iArgIdx"> (Optional) which solution to return (only 1 is implemented).</param>
	// <returns> An InferenceTermination status.</returns>
	InferenceTermination arg(std::vector<LabelType>& vLabels, const size_t iArgIdx = 1) const; 
	
	ValueType value() const; 
	void reset();
	ValueType bound() const { return dBound; };

private:
	
	GM& gmFine;  ///< The original model
	MSGM msgmFine; ///< Multiscale model of the finest level
	SubParam algParam;	///< Parameters for the relaxation algorithm
	SubParam algParamLast;	///< Parameters for the last relaxations
		
	ValueType dBound;

	Parameter param;
	double dMaxEntropy;	///< The maximum entropy (uniform on the no. of labels)

#ifdef DEBUG
	int nDepth;		///< The depth the algorithm reached
#endif

	typedef TRWSi<GM, Minimizer> TRWS_i;
	void initParam(typename TRWS_i::Parameter& algParam, int nMaxItr);
	void inline initRuntimeParam(typename TRWS_i::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef Bruteforce<GM,Minimizer> Brute;
	void initParam(typename Brute::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename Brute::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef LazyFlipper<GM, Minimizer> LF;
	void initParam(typename LF::Parameter& param, int nMaxItr);
	void initRuntimeParam(typename LF::Parameter& param, vector<LabelType> vLabels);

#ifndef _MSC_VER
	typedef FastPD<GM> FPD;
	void initParam(typename FPD::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename FPD::Parameter& param, 
		vector<LabelType> vLabels) {};
	/*typedef MRFLIB<GM> MRFLB;
	void initParam(typename MRFLB::Parameter& param);*/
    typedef DDDualVariableBlock2<marray::View<ValueType,false>> DualBlockType;
	typedef DualDecompositionBundle<GM, DynamicProgramming<
		typename DualDecompositionBase<GM, DualBlockType>::SubGmType, Minimizer>, 
		DualBlockType> DDBundle;
	void initParam(typename DDBundle::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename DDBundle::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef TRWS<GM> extTRWS;
	void initParam(typename extTRWS::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename extTRWS::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef GCOLIB<GM> GCO;
	void initParam(typename GCO::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename GCO::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef DualDecompositionSubGradient<GM, DynamicProgramming<
		typename DualDecompositionBase<GM, DualBlockType>::SubGmType, Minimizer>, 
		DualBlockType> DDSubGrad;
	void initParam(typename DDSubGrad::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename DDSubGrad::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef Multicut<GM, Minimizer> MCA;
	void initParam(typename MCA::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename MCA::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef LPCplex<GM,Minimizer> CPLEX;
	void initParam(typename CPLEX::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename CPLEX::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef InfAndFlip<GM,Minimizer,extTRWS> TRWS_LF2;
	void initParam(typename TRWS_LF2::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename TRWS_LF2::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef AlphaBetaSwap<GM,QPBO<GM>> SWAP_QPBO;
	void initParam(typename SWAP_QPBO::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename SWAP_QPBO::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef AStar<GM, Minimizer> ASTAR;
	void initParam(typename ASTAR::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename ASTAR::Parameter& param, 
		vector<LabelType> vLabels) {};
	typedef QPBO<GM> extQPBO;
	void initParam(typename extQPBO::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename extQPBO::Parameter& param, 
		vector<LabelType> vLabels);
	typedef MQPBO<GM, Minimizer> mQPBO;
	void initParam(typename mQPBO::Parameter& param, int nMaxItr);
	void inline initRuntimeParam(typename mQPBO::Parameter& param,
		vector<LabelType> vLabels);

	typedef typename SWAP_QPBO::Parameter QpboParam;
	QpboParam qpboParam;
#endif

	typedef struct EdgeScore_ {
		Edge* pEdge;
		IndexType idx;
		double dScore;
		EdgeScore_ (Edge* pEdge, IndexType idx, double dScore) :
			pEdge(pEdge), idx(idx), dScore(dScore) {}
	} EdgeScore;
#ifdef _MSC_VER
	void vCycle(bool bFinest, MSGM& msGM,
		vector<EdgeScore> vEdgeScores =
		vector<EdgeScore>());
#else
	void vCycle(bool bFinest, MSGM& msGM, 
		vector<typename MSMRF<GM,ALG>::EdgeScore> vEdgeScores = 
		vector<typename MSMRF<GM,ALG>::EdgeScore>());
#endif
	template<class ALG_> ALG_* relax(GM& gm, MSGM& msGM, 
		typename ALG_::Parameter& param, ALG_* pAlg = NULL);
	void reParameterize(MSGM& msGM);
	void scoreEdges(MSGM& msGM, vector<EdgeScore>& vEdgeScores);
	inline void shuffleEdges(vector<EdgeScore>& vEdgeScores);
	void coarsenGraph(MSGM& msGMfine, vector<EdgeScore>& vEdgeScores, 
		MSGM& msGMcoarse);
	void prolong(MSGM& msGMcoarse, MSGM& msGMfine);
	void makeCRgraph(MSGM& msGMcoarse, MSGM& msGM, MSGM& msGMcr, vector<IndexType>& vidxMapCR2fine);
	void prolongCR(MSGM& msGMcr, MSGM& msGMfine, vector<IndexType>& vidxMapCR2fine);
	void solve(MSGM& msGMcoarse);
	
	static bool compareScoreEdge(EdgeScore lhs, EdgeScore rhs) {
		return lhs.dScore < rhs.dScore;
	};

#ifdef DEBUG
	struct Unary {
		double d0;
		double d1;
	};
	struct Pairwise {
		int v0;
		int v1;
		double d00;
		double d01;
		double d10;
		double d11;
	};

	static void printModel(MSGM& msGM)
	{
		vector<Unary> vUnary;
		vUnary.resize(msGM.numberOfVariables());
		vector<Pairwise> vPairwise;
		vPairwise.resize(msGM.numberOfFactors() - vUnary.size());
		
		int iPairwiseInd = 0;
		for (int f = 0; f < msGM.numberOfFactors(); f++) {

			const typename MSGM::FactorType& fac = msGM[f];

			if (fac.numberOfVariables() == 1) {
				vUnary[fac.variableIndex(0)].d0 = fac(0);
				vUnary[fac.variableIndex(0)].d1 = fac(1);
			}
			else {
				vPairwise[iPairwiseInd].v0 = fac.variableIndex(0);
				vPairwise[iPairwiseInd].v1 = fac.variableIndex(1);

				IndexType vars[] = {0,0};
				vPairwise[iPairwiseInd].d00 = fac(vars);
				vars[1] = 1;
				vPairwise[iPairwiseInd].d01 = fac(vars);
				vars[0] = 1;
				vPairwise[iPairwiseInd].d11 = fac(vars);
				vars[1] = 0;
				vPairwise[iPairwiseInd].d10 = fac(vars);
								
				iPairwiseInd++;
			}

		}

		int i = 0;
	};
#endif

};

template<class GM, class ALG>
MSMRF<GM,ALG>::MSMRF(GM& gm, const Parameter& param) : 
	gmFine(gm),		// save a reference to original model
	msgmFine(gm),	// construct the multiscale graphical layer of the finest scale
	param(param),	// algorithm parameters
	dMaxEntropy(log((double)msgmFine.nLabels)),	// calculate max entropy according to # labels
	dBound(0)
{	
#ifdef DEBUG
	nDepth = 0;
#endif

	initParam(algParam, param.nSubInfMaxIter);
	initParam(algParamLast, param.nSubInfFinalMaxIter);
	
#ifndef _MSC_VER
	qpboParam.maxNumberOfIterations_ = param.nSubInfMaxIter;
#endif

}


template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename TRWS_i::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.maxNumberOfIterations_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM, ALG>::initParam(typename Brute::Parameter& algParam, int nMaxItr)
{
}

template<class GM, class ALG>
void MSMRF<GM, ALG>::initParam(typename LF::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
}
template<class GM, class ALG>
void MSMRF<GM, ALG>::initRuntimeParam(typename LF::Parameter& param, vector<LabelType> vLabels)
{	
}

#ifndef _MSC_VER

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename FPD::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.numberOfIterations_ = nMaxItr;
}
/*
template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename MRFLB::Parameter& algParam)
{
	algParam = param.subInfParam;
	algParam.numberOfIterations_ = param.nSubInfMaxIter;
	algParam.inferenceType_ = MRFLB::Parameter::EXPANSION;
	algParam.energyType_ = MRFLB::Parameter::TABLES;
}
*/

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename DDBundle::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.maximalNumberOfIterations_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename extTRWS::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.numberOfIterations_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename GCO::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.numberOfIterations_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename DDSubGrad::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.maximalNumberOfIterations_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename MCA::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.timeOut_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename CPLEX::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.timeLimit_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename TRWS_LF2::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.subPara_.numberOfIterations_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename SWAP_QPBO::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.maxNumberOfIterations_ = nMaxItr;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename ASTAR::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	// TODO: AStar does not have a max-iter property only objective bound which is an
	// energy value to stop at
	// algParam.maxNumberOfIterations_ = nMaxItr; 
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::initParam(typename extQPBO::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	algParam.useImproveing_ = true;
}
template<class GM, class ALG>
void inline MSMRF<GM,ALG>::initRuntimeParam(typename extQPBO::Parameter& param, 
											vector<LabelType> vLabels)
{
	param.label_ = vLabels;
}


template<class GM, class ALG>
void MSMRF<GM, ALG>::initParam(typename mQPBO::Parameter& algParam, int nMaxItr)
{
	algParam = param.subInfParam;
	//algParam.useImproveing_ = true;
}
template<class GM, class ALG>
void inline MSMRF<GM, ALG>::initRuntimeParam(typename mQPBO::Parameter& param,
	vector<LabelType> vLabels)
{
	param.label_ = vLabels;
}

#endif


template<class GM, class ALG>
MSMRF<GM,ALG>::~MSMRF(void)
{
}

// to infer the solution we invoke the multiscale recursion on the finest scale
template<class GM, class ALG>
template<class VisitorType> 
InferenceTermination MSMRF<GM,ALG>::infer(VisitorType& visitor)
{
	visitor.begin(*this);

	if (param.bImprove) {
		if (param.sImproveAlg == "BPS"){
			typedef TRWS<GM> ImproveAlg;
			typename ImproveAlg::Parameter improveAlgParam;
			improveAlgParam.doBPS_ = true;
			improveAlgParam.energyType_ = ImproveAlg::Parameter::EnergyType::TABLES;
			ImproveAlg improveAlg(msgmFine, improveAlgParam);
			improveAlg.infer();
			improveAlg.arg(msgmFine.vLabels);
			msgmFine.bLabels = true;
		}
		else if (param.sImproveAlg == "LBP") {
			typedef BeliefPropagationUpdateRules<GM, Minimizer> UpdateRulesType;
			typedef MessagePassing<GM, Minimizer, UpdateRulesType> ImproveAlg;
			typename ImproveAlg::Parameter improveAlgParam;
			improveAlgParam.maximumNumberOfSteps_ = 1000;
			improveAlgParam.damping_ = 0.5;
			improveAlgParam.bound_ = 0.00001;
			ImproveAlg improveAlg(msgmFine, improveAlgParam);
			improveAlg.infer();
			improveAlg.arg(msgmFine.vLabels);
			msgmFine.bLabels = true;
		}
		//else if (param.sImproveAlg == "EXPAND") {
		//	typedef GCOLIB<GM> ImproveAlg;
		//	typename ImproveAlg::Parameter improveAlgParam;
		//	improveAlgParam.inferenceType_ = SubInf::Parameter::InferenceType::EXPANSION;
		//	improveAlgParam.energyType_ = SubInf::Parameter::EnergyType::VIEW;
		//}
		else if (param.sImproveAlg == "ILP") {
			typedef LPCplex<GM, Minimizer> ImproveAlg;
			typename ImproveAlg::Parameter algParam;
			algParam.numberOfThreads_ = 1;
			algParam.integerConstraint_ = true;
		}
	}


	/** reparametrize **/
	reParameterize(msgmFine);
	
#ifdef DEBUG
	/*msgmFine.vLabels.resize(msgmFine.numberOfVariables());
	msgmFine.bLabels = true;
	double dEnrgy = msgmFine.evaluate(msgmFine.vLabels);*/
	//printModel(msgmFine);
#endif

#ifdef DEBUG
	/*msgmFine.vLabels.resize(msgmFine.numberOfVariables());
	msgmFine.bLabels = true;
	double dEnrgy2 = msgmFine.evaluate(msgmFine.vLabels);*/
	//printModel(msgmFine);
#endif

	/** score edges **/
	// the list of energies is a vector of tuples of directed edges and energy values.
	// the direction is identified such that the integer (2nd value in the tuple) is the
	// index of the interpolator. 
	vector<typename MSMRF<GM,ALG>::EdgeScore> vEdgeScores;
	scoreEdges(msgmFine, vEdgeScores);

	//call the recursive function on the finest scale
	for (int i = 0; i < param.nVcycles; i++) {
		
		vCycle(true, msgmFine, vEdgeScores);

		cout << "<< | V-cycle: " << setw(5) << (i+1)
			<< " | Energy: " << setw(8) << msgmFine.evaluate(msgmFine.vLabels) << endl;
		cout << "---------------------" << endl;
		
		visitor(*this);
	}

	//final relaxations
	if (this->param.nSubInfFinalMaxIter > 0)
		delete relax<ALG>(gmFine, msgmFine, algParamLast);

	cout << "<< | Final Energy: " << setw(8) << msgmFine.evaluate(msgmFine.vLabels) << endl;
	
	visitor.end(*this);

	return NORMAL;
}

template<class GM, class ALG>
InferenceTermination MSMRF<GM,ALG>::arg(vector<LabelType>& vStates, const size_t iArgIdx) const
{
	vStates = msgmFine.vLabels;
	return NORMAL;
}

template<class GM, class ALG>
typename MSMRF<GM,ALG>::ValueType MSMRF<GM,ALG>::value() const 
{
	return msgmFine.evaluate(msgmFine.vLabels.begin());
} 
 

   
// <summary> The multiscale recursive function</summary>
// <remarks> This function is called at each scale on a multiscale graphical layer
// 			 </remarks>
// <param name="bFinest">  Whether the call is made at the finest scale</typeparam>
// <param name="msGM"> [in,out] The multiscale graphical layer to operate on.</param> 
template<class GM, class ALG>
void MSMRF<GM,ALG>::vCycle(bool bFinest, MSGM& msGM, 
						   vector<EdgeScore> vEdgeScores)
{
	/** solve **/

	// check whether we reached the coarsest level
	if (msGM.vVertices.size() <= param.nSolveSize) {
#ifdef DEBUG
	cout << "** | Solved " << endl;
#endif	
		//solve
		solve(msGM);
		return;
	} 
	
	/** relaxation **/
	Inference<GM, Minimizer>* pAlg = NULL;

	// infer using the sub algorithm (only after first solve)
	if (msGM.bLabels && this->param.nSubInfMaxIter > 0) {
		if (!bFinest && param.midScaleAlg == Parameter::MID_SCALE_ALG::SWAP) {
			pAlg = relax<SWAP_QPBO>(msGM, msGM, qpboParam);
			delete pAlg;
		}
		
		else if (bFinest || param.midScaleAlg == Parameter::MID_SCALE_ALG::SAME)
			pAlg = relax<ALG>(bFinest ? gmFine : msGM, msGM, algParam);
	}

	//#ifdef DEBUG
	cout << "." << endl;
	cout << "<< Before Reparam << | Vertices: " << setw(5) << msGM.vVertices.size()
		<< " | Edges: " << setw(8) << msGM.vEdges.size();
	if (msGM.bLabels)
		cout << " | Energy: " << setw(8) << msGM.evaluate(msGM.vLabels) << endl;
	else
		cout << endl;
	//#endif

#ifdef DEBUG
	/*msGM.vLabels.resize(msGM.numberOfVariables());
	msGM.bLabels = true;
	double dEnrgy = msGM.evaluate(msGM.vLabels);*/
	//printModel(msGM);
#endif

	if (vEdgeScores.size() == 0) {

		/** reparametrize **/
		reParameterize(msGM);

#ifdef DEBUG
		/*msGM.vLabels.resize(msGM.numberOfVariables());
		msGM.bLabels = true;
		double dEnrgy2 = msGM.evaluate(msGM.vLabels);*/
		//printModel(msGM);
#endif

		/** score edges **/
		// the list of energies is a vector of tuples of directed edges and energy values.
		// the direction is identified such that the integer (2nd value in the tuple) is the
		// index of the interpolator. 
		scoreEdges(msGM, vEdgeScores);

	}

	#ifdef DEBUG
	cout << "<< After Reparam  << | Vertices: " << setw(5) << msGM.vVertices.size()
		<< " | Edges: " << setw(8) << msGM.vEdges.size();
	if (msGM.bLabels)
		cout << " | Energy: " << setw(8) << msGM.evaluate(msGM.vLabels) << endl;
	else
		cout << endl;
	#endif
	    
	// add randomization when eanbled and after initial labeling
	if (param.nBins > 0 && msGM.bLabels) {
		shuffleEdges(vEdgeScores);
		sort(vEdgeScores.begin(), vEdgeScores.end(), &compareScoreEdge);
	}

	/** coarsen graph **/  
	MSGM msGMcoarse(msGM.nLabels);
	coarsenGraph(msGM, vEdgeScores, msGMcoarse);

	//#ifdef DEBUG
	cout << "<< After Coarsen  << | Vertices: " << setw(5) << msGMcoarse.vVertices.size()
		<< " | Edges: " << setw(8) << msGMcoarse.vEdges.size();
	if (msGM.bLabels)
		cout << " | Energy: " << setw(8) << msGMcoarse.evaluate(msGMcoarse.vLabels) << endl;
	else
		cout << endl;
	//#endif

#ifdef DEBUG
	//SWAP_QPBO alg2(msGMcoarse);
	//relax(msGMcoarse, alg2);
	//double ddd = msGMcoarse.evaluate(msGMcoarse.vLabels.begin());

	//if (msGMcoarse.vVertices.size() <= 2)
	//	printModel(msGMcoarse);
	/*vector<double> vlbls(msGMcoarse.vVertices.size());
	double dCoarseEnrgy = msGMcoarse.evaluate(vlbls);
	vlbls.resize(msGM.vVertices.size());
	double dFineEnrgy = msGM.evaluate(vlbls);
	assert(abs(dCoarseEnrgy - dFineEnrgy) < 0.0001);*/
#endif

	/** recursive call **/
	vCycle(false, msGMcoarse);
	
	
	#ifdef DEBUG
	cout << "<< Before Prolong << | Vertices: " << setw(5) << msGMcoarse.vVertices.size()
		<< " | Edges: " << setw(8) << msGMcoarse.vEdges.size()
		<< " | Energy: " << setw(8) << msGMcoarse.evaluate(msGMcoarse.vLabels) << endl;
	#endif

	/** prolongation **/
	if (param.bCmptRelax) {
		// doing compatible relaxation

		// make CR graph
		MSGM msGMcr(msGM.nLabels);
		vector<IndexType> vidxMapCR2fine;
		makeCRgraph(msGMcoarse, msGM, msGMcr, vidxMapCR2fine);

		// relax CR graph
		pAlg = relax<ALG>(msGMcr, msGMcr, algParam, (ALG*)pAlg);
		delete pAlg;

		// prolong CR labels
		// prolongCR(msGMcr, msGM, vidxMapCR2fine);
		// prolong(msGMcoarse, msGM);

	}
	else {
		// usual prolongation
		prolong(msGMcoarse, msGM);
	}
	

	#ifdef DEBUG
	cout << "<< After Prolong  << | Vertices: " << setw(5) << msGM.vVertices.size()
		<< " | Edges: " << setw(8) << msGM.vEdges.size()
		<< " | Energy: " << setw(8) << msGM.evaluate(msGM.vLabels) << endl;

	#endif

#ifdef DEBUG
	//printModel(msGM);
	double d1 = msGMcoarse.evaluate(msGMcoarse.vLabels);
	double d2 = msGM.evaluate(msGM.vLabels);
	double d = d2 - d1;
	assert(d < 1e-6);
#endif

	/** relaxation **/
	// infer using the sub algorithm 
	if (this->param.nSubInfMaxIter > 0) {
		if (!bFinest && param.midScaleAlg == Parameter::MID_SCALE_ALG::SWAP) 
			pAlg = relax<SWAP_QPBO>(msGM, msGM, qpboParam);

		else if (bFinest || param.midScaleAlg == Parameter::MID_SCALE_ALG::SAME)
			//pAlg = relax<ALG>(bFinest ? gmFine : msGM, msGM, algParam, (ALG*)pAlg);
			pAlg = relax<ALG>(bFinest ? gmFine : msGM, msGM, algParam, (ALG*)pAlg);
	}
	if (pAlg)
		delete pAlg;

	#ifdef DEBUG
	cout << "<<  After Relax   << | Vertices: " << setw(5) << msGM.vVertices.size()
		<< " | Edges: " << setw(8) << msGM.vEdges.size()
		<< " | Energy: " << setw(8) << msGM.evaluate(msGM.vLabels) << endl
		<< "." << endl;
	#endif
}

// <summary> Executes the single scale algorithm on the current model.</summary>
// <remarks> The algorithm is assumed to be initialized with the proper parameters.</remarks>
// <param name="msGM"> [in,out] The multiscale graphical layer to operate on.</param>
// <param name="alg">  [in,out] The algorithm instance to use for relaxation.</param>
template<class GM, class ALG>
template<class ALG_>
ALG_* MSMRF<GM,ALG>::relax(GM& gm, MSGM& msGM, typename ALG_::Parameter& param,
						  ALG_* pAlg)
{
#ifdef DEBUG

	//for (int v = 0; v < msGM.vLabels.size(); v++)
	//	msGM.vLabels[v] = floor((double)rand() / (double)RAND_MAX + 0.5);

#endif
	
	// some solvers (like QPBO) require the labeling in param
	if (!pAlg) {
		initRuntimeParam(param, msGM.vLabels);

		// instantiate the solver
		pAlg = new ALG_(gm, param);
	}

	// set the alg's starting point to the current states
	pAlg->setStartingPoint(msGM.vLabels.begin());

	// infer using the algorithm
	pAlg->infer();

	// update the current model with the argmin
	pAlg->arg(msGM.vLabels);

	// iterate all vertices to update post-predictive value
	for (int i = 0; i < msGM.vVertices.size(); i++) {
		Vertex*& pv = msGM.vVertices[i];

		// increase the matching entry by 1
		(*pv->pvdPostPred)[msGM.vLabels[pv->idx]]++;
		(*pv->pdPostPredSum)++;
	}

	// touched the labels
	msGM.bLabels = true;

	dBound = pAlg->bound();

	return pAlg;
}


// <summary> Re parameterize.</summary>
// <remarks> Stav, 06/01/2014.</remarks>
// <param name="msGM"> [in,out] The milliseconds gm.</param>
template<class GM, class ALG>
void MSMRF<GM,ALG>::reParameterize(MSGM& msGM)
{
	double dVal;			//current function value
	LabelType lbl[2];		//label indexers

#ifdef DEBUG
	int n = msGM.vEdges.size();
#endif

#ifdef DEBUG
	/*msGM.vLabels.resize(msGM.numberOfVariables());
	msGM.bLabels = true;
	double dEnrgy = msGM.evaluate(msGM.vLabels);*/
	int nFuncs = 0;
	set<FunctionIdentifier> setFn;
#endif
	
	//marginal means
	struct Means {
		vector<double> vdMeanLeft;
		vector<double> vdMeanRight;

		Means() {}
		Means(int nLabels) :
			vdMeanLeft(nLabels),
			vdMeanRight(nLabels)
		{}
	};

	typedef map<FunctionIdentifier, Means> MeansMap;
	typedef typename MeansMap::iterator ItrMeans;
	MeansMap mfivdMeans;
	
	// iterate all edges (dual var factors)
	for (int i = 0; i < msGM.vEdges.size(); i++) {
		Edge& edge = *msGM.vEdges[i];	//the current edge in question
		Means* pMeans;					//the means of the function behind this edge

		// see if the function behind this factor was already updated
		ItrMeans itrMeans = mfivdMeans.find(edge.fi);
		if (itrMeans == mfivdMeans.end()) {
			Means means(msGM.nLabels);
			
			// calculate marginal means
			//fill(means.vdMeanLeft.begin(), vdMeanLeft.end(), 0);
			//fill(vdMeanRight.begin(), vdMeanRight.end(), 0);

			for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {
				for (lbl[1] = 0; lbl[1] < msGM.nLabels; lbl[1]++) {
					dVal = edge(lbl);
					means.vdMeanLeft[lbl[0]] += (dVal / msGM.nLabels);
					means.vdMeanRight[lbl[1]] += (dVal / msGM.nLabels);
#ifdef DEBUG
					double dfn = edge.fn()(lbl[1], lbl[0]);
					double dfac = edge(lbl);
					assert(edge(lbl) == edge.fn()(lbl[0], lbl[1]));
#endif
				}
			}

			//update pairwise function
			for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {
				for (lbl[1] = 0; lbl[1] < msGM.nLabels; lbl[1]++) {
					edge.fn()(lbl[0], lbl[1]) -= 
						(means.vdMeanLeft[lbl[0]] + means.vdMeanRight[lbl[1]]);
				}
			}

			mfivdMeans[edge.fi] = means;
			pMeans = &mfivdMeans[edge.fi];

#ifdef DEBUG
			nFuncs++;
#endif
		}

		// if the function was already reparameterized get the updated means
		else
			pMeans = &(itrMeans->second);
		
		//update the unaries
		for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {
			// update unary term for left and right (we use the left indexer 
			// but this has no meaning) 
			edge.v1.fn()(lbl[0]) += pMeans->vdMeanLeft[lbl[0]];
			edge.v2.fn()(lbl[0]) += pMeans->vdMeanRight[lbl[0]];
		}

#ifdef DEBUG
		//assert(setFn.find(edge.v1.fi) == setFn.end());
		//assert(setFn.find(edge.v2.fi) == setFn.end());
		setFn.insert(edge.v1.fi);
		setFn.insert(edge.v2.fi);
#endif

#ifdef DEBUG
		//double dEnrgyEdge = msGM.evaluate(msGM.vLabels);
		//assert(abs(dEnrgyEdge - dEnrgy) < 0.0001);
#endif

	}	// factors (edges) loop
	
}

// <summary> Scores the edges (dual-var factors) for coarsening selection.</summary>
// <remarks> This function iterates all the edges and their various assignments. The edges are
// treated in a directed manner so that the list scores both directions for each edge.</remarks>
// <param name="msGM">		  [in,out] The multiscale graphical model to operate on.</param>
// <param name="vEdgeScores"> [in,out] The edge scores
template<class GM, class ALG>
void MSMRF<GM,ALG>::scoreEdges(MultSclGM<GM>& msGM, 
							   vector<EdgeScore>& vEdgeScores)
{
	//reserve space for the scores (times 2 for both directions)
	vEdgeScores.reserve(msGM.numberOfFactors() * 2);

	vector<vector<double>> vvdNormUnary(msGM.numberOfVariables(), 
		vector<double>(msGM.nLabels));	//nornalized unary terms
	
	double dVal;							//current, mean and sum of all values of edge
	double dUnaryLeft, dUnaryRight;			//unary terms of left and right variables
	double dSumPr;							//current probabilities sums, for normalizing
	
	vector<double> vdPrLeft(msGM.nLabels);		//(unnormalized) margninal probabilities
	vector<double> vdPrRight(msGM.nLabels);
	vector<double> vdEntLeft(msGM.nLabels);	//entropy
	vector<double> vdEntRight(msGM.nLabels);

	double dHLeft, dHRight;					//intermidates (see below)
	LabelType lbl[2];						//a possible assignment 0- left, 1- right
	
	// iterate and normalize unary terms
	if (param.bDegEnt) {
		for (int i = 0; i < msGM.vVertices.size(); i++) {
			Vertex*& pv = msGM.vVertices[i];

			// iterate labels
			for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {
				// minus 1 on the number of factors is the degree (excluding the unary term)
				vvdNormUnary[pv->idx][lbl[0]] = (*pv)(lbl[0]) / 
					(msGM.numberOfFactors(pv->idx) - 1);
			}
		}
	}

	// iterate all edges (dual var factors)
	for (int i = 0; i < msGM.vEdges.size(); i++) {
		Edge& edge = *msGM.vEdges[i];

		dSumPr = 0;

		//zero out the marginal probabilities and entropies
		fill(vdPrLeft.begin(), vdPrLeft.end(), 0);
		fill(vdPrRight.begin(), vdPrRight.end(), 0);
		fill(vdEntLeft.begin(), vdEntLeft.end(), 0);
		fill(vdEntRight.begin(), vdEntRight.end(), 0);
		
		// first pass - find average of fused binary and unary terms.
		// iterate labels of the left vertex.
		for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {
			
			//iterate labels of right vertex
			for (lbl[1] = 0; lbl[1] < msGM.nLabels; lbl[1]++) {

				dVal = edge(lbl);

				// get unary terms (either direct or normalized)
				if (param.bDegEnt) {
					dUnaryLeft = vvdNormUnary[edge.v1.idx][lbl[0]];
					dUnaryRight = vvdNormUnary[edge.v2.idx][lbl[1]];
				}
				else {
					dUnaryLeft = edge.v1(lbl[0]);
					dUnaryRight = edge.v2(lbl[1]);
				}

				// add the corresponding unary terms of both variables
				// and map to probability space
				dVal += dUnaryLeft + dUnaryRight;
				dVal = exp(-dVal);
			
				// update probabilities and entropy
				dSumPr += dVal;
				vdPrLeft[lbl[0]] += dVal;
				vdPrRight[lbl[1]] += dVal;
				vdEntLeft[lbl[0]] += -dVal * log(dVal);
				vdEntRight[lbl[1]] += -dVal * log(dVal);
				
			}
		}
		
		// update the enropy to account for the normalization factor
		dHLeft = 0;
		dHRight = 0;
		
		// check whether working with posterior predictive
		if (param.bUsePostPred) {
			for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {

				// update entropy for left vertex
				 dHLeft += (*edge.v1.pvdPostPred)[lbl[0]] / (*edge.v1.pdPostPredSum) *
					(vdEntLeft[lbl[0]] / vdPrLeft[lbl[0]] + log(vdPrLeft[lbl[0]]));
				
				// update entropy for right vertex
				dHRight += (*edge.v2.pvdPostPred)[lbl[0]] / (*edge.v2.pdPostPredSum) *
					(vdEntRight[lbl[0]] / vdPrRight[lbl[0]] + log(vdPrRight[lbl[0]]));
			}
		}
		// not working with the post. predictive
		else {
			for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {

				// update entropy for left vertex
				dHLeft += (vdEntLeft[lbl[0]] + vdPrLeft[lbl[0]] * 
					log(vdPrLeft[lbl[0]])) / dSumPr;
			
				// update entropy for right vertex
				dHRight += (vdEntRight[lbl[0]] + vdPrRight[lbl[0]] * 
					log(vdPrRight[lbl[0]])) / dSumPr;
			}
		}
		
		// add the score to the list
		vEdgeScores.emplace_back(&edge, 0, dHLeft);
		vEdgeScores.emplace_back(&edge, 1, dHRight);
		
	}

	// sort the list according to entropy
	sort(vEdgeScores.begin(), vEdgeScores.end(), &compareScoreEdge);

	// bin edges
	if (param.nBins > 0) {
		for (int i = 0; i < vEdgeScores.size(); i++) {
			// bin edges - normalize by the maximal entropy
			vEdgeScores[i].dScore /= dMaxEntropy;
			
			// round to nearest bin
			vEdgeScores[i].dScore = floor(vEdgeScores[i].dScore * param.nBins + 0.5);
		}
	}
}

template<class GM, class ALG>
inline void MSMRF<GM,ALG>::shuffleEdges(vector<EdgeScore>& vEdgeScores)
{
	// iterate scores and add randomness for intra-bin ordering
	for (int i = 0; i < vEdgeScores.size(); i++)
		vEdgeScores[i].dScore += 0.5 * (double)rand() / (double)RAND_MAX;
}

template<class GM, class ALG>
void MSMRF<GM,ALG>::coarsenGraph(MSGM& msGMfine, vector<EdgeScore>& vEdgeScores, 
								 MSGM& msGMcoarse)
{
	IndexType shp[2] = {msGMcoarse.nLabels, msGMcoarse.nLabels};
	LabelType lbl[2];		//a possible assignment 0- left, 1- right
	int nSkipped = 0;
	double d;
	
	//reset the user variable for all edges
	for_each(msGMfine.vEdges.begin(), msGMfine.vEdges.end(), 
		[](Edge*& pEdge) { pEdge->dwUser = FALSE; });

	// allocate some space for the coarse model. 
	msGMcoarse.reserveVertices((int)( (float)
		msGMfine.vVertices.size() * COARSE_VERTEX_RATIO ));

	// first pass - iterate over edges according to score
	for (int i = 0; i < vEdgeScores.size(); i++) {

		EdgeScore& edgeScore = vEdgeScores[i];
		Edge* pEdge = edgeScore.pEdge;
		IndexType iLeft = edgeScore.idx;
	
		//check edge has not been processed (in any direction)
		if (!pEdge->dwUser) {

			// get pointers according to the selected edge direction
			bool bRev;				//"reversed" edge?
			Vertex *pvLeft, *pvRight;
			if (iLeft == 0) {
				pvLeft = &pEdge->v1;
				pvRight = &pEdge->v2;
				bRev = false;
			}
			else {
				pvLeft = &pEdge->v2;
				pvRight = &pEdge->v1;
				bRev = true;
			}
		
		
			// is the right vertex still available (is not part of the coarse graph yet)
			if (!pvRight->pvCoarse) {

				// is the left vertex still available
				if (!pvLeft->pvCoarse) {
				
					ExFunction fnUnary(shp, shp+1);
					double dMinEnergy, dEnergy;
					LabelType lblArgmin;

					//create a new coarse vertex		
					Vertex& vCoarse = msGMcoarse.addVertex(pvLeft); 

					// maintain linking
					vCoarse.addRight(pvRight);

					// iterate labels of the left vertex
					for (lbl[0] = 0; lbl[0] < msGMfine.nLabels; lbl[0]++) {
						
						// if this is the currently assigned label to the fine vertex
						if (msGMfine.bLabels && msGMfine.vLabels[pvLeft->idx] == lbl[0]) {
							// in this case, we maintain monotonicity by making sure
							// prolongation will recover the current assignment
							lblArgmin = msGMfine.vLabels[pvRight->idx];
							lbl[1] = lblArgmin;
							dMinEnergy = (*pvRight)(lblArgmin) + (*pEdge)(lbl,bRev);

							// also, let current assignment in the coarse graph be equal
							// we resize before to make sure we do not overflow
							msGMcoarse.vLabels.resize(vCoarse.idx + 1);
							msGMcoarse.vLabels[vCoarse.idx] = lbl[0];
						}

						// no assigned labels
						else {
							dMinEnergy = numeric_limits<double>::infinity();
					
							// find the minimum assignment for the right vertex
							for (lbl[1] = 0; lbl[1] < msGMfine.nLabels; lbl[1]++) {
								dEnergy = (*pEdge)(lbl,bRev) + (*pvRight)(lbl[1]);

								if (dEnergy < dMinEnergy) {
									dMinEnergy = dEnergy;
									lblArgmin = lbl[1];
								}
							}
						}

						// set the entry for the right vertex in the prolongation table
						pvRight->vlblPlongTab[lbl[0]] = lblArgmin;

						// set the unary value for the left vertex
						fnUnary(lbl[0]) = dMinEnergy + (*pvLeft)(lbl[0]);
					}

					// set the corase vertex unary factor
					vCoarse.setUnary(fnUnary);
				
					pEdge->dwUser = TRUE;

				}

				// is the left vertex in the coarse graph but as an interpolator
				else if (pvLeft->isLeft()) {

					double dMinEnergy, dEnergy;
					LabelType lblArgmin;
					Vertex* pvCoarse = pvLeft->pvCoarse;

					// maintain linking (add another right vertetx)
					pvCoarse->addRight(pvRight);

					// iterate labels of the left vertex
					for (lbl[0] = 0; lbl[0] < msGMfine.nLabels; lbl[0]++) {

						// if this is the currently assigned label to the fine vertex
						if (msGMfine.bLabels && msGMfine.vLabels[pvLeft->idx] == lbl[0]) {
							// in this case, we maintain monotonicity by making sure
							// prolongation will recover the current assignment
							lblArgmin = msGMfine.vLabels[pvRight->idx];
							lbl[1] = lblArgmin;
							dMinEnergy = (*pvRight)(lblArgmin) + (*pEdge)(lbl,bRev);
						}

						// no assigned labels
						else {
							dMinEnergy = numeric_limits<double>::infinity();
					
							// find the minimum assignment for the right vertex
							for (lbl[1] = 0; lbl[1] < msGMfine.nLabels; lbl[1]++) {
								dEnergy = (*pEdge)(lbl,bRev) + (*pvRight)(lbl[1]);

								if (dEnergy < dMinEnergy) {
									dMinEnergy = dEnergy;
									lblArgmin = lbl[1];
								}
							}
						}

						// set the entry for the right vertex in the prolongation table
						pvRight->vlblPlongTab[lbl[0]] = lblArgmin;

						// set the unary value for the left vertex
						d = pvCoarse->fn()(lbl[0]);
						pvCoarse->fn()(lbl[0]) += dMinEnergy;
						d = pvCoarse->fn()(lbl[0]);
					}
				
					pEdge->dwUser = TRUE;
				}
				else {
					//skipped edge
					nSkipped++;
				}

			}
			else {
				//skipped edge
				nSkipped++;
			}
		}

		// if we coarsened a labeled graph it is also labeled
		msGMcoarse.bLabels = msGMfine.bLabels;
			
	}

	// adjacency map for the corase graph 
	// TODO: find a way to preallocate enough space, we know
	// that it will be at most the number of edges in the fine graph. 
	typedef pair<Vertex*, Vertex*> VertexPair;
	typedef map<VertexPair, Edge*> AdjMap;
	AdjMap mAdj;
	IndexType pidx[2];
	
	msGMcoarse.reserveEdges(nSkipped);

	// second pass - iterate skipped edges
	for (int i = 0; i < msGMfine.vEdges.size(); i++) {
		Edge& edge = *msGMfine.vEdges[i];

		if (!edge.dwUser) {
			
			//if one of the vertices does not have a coarse vertex it is an orphan
			if (!edge.v1.pvCoarse) {
				//create a new coarse vertex		
				Vertex& vCoarse = msGMcoarse.addVertex(&edge.v1); 

				//copy the unary term
				vCoarse.setUnary(edge.v1.fn());

				if (msGMfine.bLabels) {
					// also, let current assignment in the coarse graph be equal
					// we resize before to make sure we do not overflow
					msGMcoarse.vLabels.resize(vCoarse.idx + 1);
					msGMcoarse.vLabels[vCoarse.idx] = msGMfine.vLabels[edge.v1.idx];
				}
			}
			else if (!edge.v2.pvCoarse) {
				Vertex& vCoarse = msGMcoarse.addVertex(&edge.v2); 
				vCoarse.setUnary(edge.v2.fn());
				if (msGMfine.bLabels) {
					msGMcoarse.vLabels.resize(vCoarse.idx + 1);
					msGMcoarse.vLabels[vCoarse.idx] = msGMfine.vLabels[edge.v2.idx];
				}
			}

			//check if edge is a self loop
			if (edge.v1.pvCoarse == edge.v2.pvCoarse) {
				
				//update the coarse unary term with possible assignments
				for (LabelType lblCoarse = 0; lblCoarse < msGMfine.nLabels; lblCoarse++) {
					lbl[0] = edge.v1.vlblPlongTab[lblCoarse];
					lbl[1] = edge.v2.vlblPlongTab[lblCoarse];

					//add the pairwise cost in this assignment
					//d = edge.v1.pvCoarse->fn()(lbl[0]);
					edge.v1.pvCoarse->fn()(lblCoarse) += edge(lbl);
					//d = edge.v1.pvCoarse->fn()(lbl[0]);
				}
			}
			else {

				//check if edge exists in coarse graph
				const VertexPair& vpCoarseEdge = (edge.v1.pvCoarse->idx < edge.v2.pvCoarse->idx) ? 
					make_pair(edge.v1.pvCoarse, edge.v2.pvCoarse) :
					make_pair(edge.v2.pvCoarse, edge.v1.pvCoarse);
				//const VertexPair& vpCoarseEdge = 
				//	make_pair(edge.v1.pvCoarse, edge.v2.pvCoarse);
#ifdef _MSC_VER
				typedef AdjMap::iterator ItrAdj;
#else
				typedef typename AdjMap::iterator ItrAdj;
#endif
				ItrAdj itrAdj = mAdj.find(vpCoarseEdge);
				bool bNew = (itrAdj == mAdj.end());

				// take the existing function or create a new one according to whether the edge is "new". 
				ExFunction fnTemp(shp, shp+2);
				ExFunction& fnBinary = bNew ? fnTemp : itrAdj->second->fn();
				
				// choose the correct pairwise calculation. 
				function<ValueType(LabelType*)> fCalc;
				Vertex* pvRight;
				
				// the case where both vertices are "right"
				if (!edge.v1.isLeft() && !edge.v2.isLeft()) {
					fCalc = [&](LabelType* lbl) {
						LabelType lblFine[2];
						lblFine[0] = edge.v1.vlblPlongTab[lbl[0]];
						lblFine[1] = edge.v2.vlblPlongTab[lbl[1]];
						return edge(lblFine);
					};
				}

				//the case when only one of them is right
				else if (edge.v1.isLeft() != edge.v2.isLeft()) {
					//check which one is left
					IndexType iLeft;
					if (edge.v1.isLeft()) {
						iLeft = 0;
						pvRight = &edge.v2;
					}
					else {
						iLeft = 1;
						pvRight = &edge.v1;
					}

					fCalc = [&](LabelType* lbl) {
						LabelType lblFine[2];
						lblFine[iLeft] = lbl[iLeft];
						lblFine[!iLeft] = pvRight->vlblPlongTab[lbl[!iLeft]];
						return edge(lblFine);
					};
				}

				//the case when both are left
				else {
					fCalc = [&](LabelType* lbl) {
						return edge(lbl);
					};
				}
				

				//reassign pairwise value to the coarse graph
				LabelType lbl_[2];
				for (lbl[0] = 0; lbl[0] < msGMcoarse.nLabels; lbl[0]++) {
					for (lbl[1] = 0; lbl[1] < msGMcoarse.nLabels; lbl[1]++) {
						d = fCalc(lbl);
						d = fnBinary(lbl);	
						//check pairwise orientation in coarse and fine graph
						if (vpCoarseEdge.first->idx == edge.v1.pvCoarse->idx) {
							//pairwise is oriented in the same direction
							lbl_[0] = lbl[0];
							lbl_[1] = lbl[1];
						}
						else {
							//pairwise is oriented in different directions
							lbl_[0] = lbl[1];
							lbl_[1] = lbl[0];
						}
						//update the pairwise values
						if (bNew)
							fnBinary(lbl) = fCalc(lbl_);
						else
							fnBinary(lbl) += fCalc(lbl_);
						d = fnBinary(lbl);
					}
				}

				//if this is a "new" coarse edge create it in the model
				if (bNew) {
					if (edge.v1.pvCoarse->idx < edge.v2.pvCoarse->idx) {
						pidx[0] = edge.v1.pvCoarse->idx;
						pidx[1] = edge.v2.pvCoarse->idx;
					}
					else {
						pidx[0] = edge.v2.pvCoarse->idx;
						pidx[1] = edge.v1.pvCoarse->idx;
					}
					Edge& edgeCoarse = msGMcoarse.addEdge(
						msGMcoarse.addFunction(fnBinary), pidx, pidx+2);
					mAdj[vpCoarseEdge] = &edgeCoarse;
					//mAdj[make_pair(vpCoarseEdge.second, vpCoarseEdge.first)] = &edgeCoarse;
				}

			}
		}
	}

	// initialize the labels vectors
	msGMcoarse.vLabels.resize(msGMcoarse.vVertices.size());
}

//Solve - direct solution to a very small graph
template<class GM, class ALG>
void MSMRF<GM,ALG>::solve(MSGM& msGM)
{
	// in case of 1 variable
	if (msGM.vVertices.size() == 1) {
		// choose minimum energy label
		double dMinEnergy = numeric_limits<double>::infinity();

		for (LabelType lbl = 0; lbl < msGM.nLabels; lbl++) {
			double dEnergy = (*msGM.vVertices[0])(lbl);
			if (dEnergy < dMinEnergy) {
				msGM.vLabels[0] = lbl;
				dMinEnergy = dEnergy;
			}
		}
	}

	// 2 vertices remaining
	else if (msGM.vVertices.size() == 2) {
		//find the minium pairwise assignment
		double dMinEnergy = numeric_limits<double>::infinity();
		LabelType lbl[2];
		Edge& edge = *msGM.vEdges.front();

		for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++) {
			for (lbl[1] = 0; lbl[1] < msGM.nLabels; lbl[1]++) {
				
				double dEnergy = edge.v1(lbl[0]) + edge.v2(lbl[1]) + edge(lbl);

				if (dEnergy < dMinEnergy) {
					msGM.vLabels[edge.v1.idx] = lbl[0];
					msGM.vLabels[edge.v2.idx] = lbl[1];
					dMinEnergy = dEnergy;
				}
			}
		}

	}

	// many vertices - only if labels are not initialized
	else if (!msGM.bLabels) {

		for (int v = 0; v < msGM.vVertices.size(); v++) {
			double dMinEnergy = numeric_limits<double>::infinity();

			for (LabelType lbl = 0; lbl < msGM.nLabels; lbl++) {
				double dEnergy = (*msGM.vVertices[v])(lbl);
				if (dEnergy < dMinEnergy) {
					msGM.vLabels[v] = lbl;
					dMinEnergy = dEnergy;
				}
			}
		}

		//Bruteforce<GM, Minimizer> brute(msGM);
		//brute.infer();
		//brute.arg(msGM.vLabels,0);

		//LPCplex<GM, Minimizer> ilp(msGM);
		//ilp.infer();
		//ilp.arg(msGM.vLabels);
	}

	// we touched the labels
	msGM.bLabels = true;
}


//Prolong
template<class GM, class ALG>
void MSMRF<GM,ALG>::prolong(MSGM& msGMcoarse, MSGM& msGMfine)
{
	
	for (int i = 0; i < msGMfine.vVertices.size(); i++) {
		Vertex*& pv = msGMfine.vVertices[i];

		// if we are the left vertex in the coarse graph
		if (pv->isLeft())
			msGMfine.vLabels[pv->idx] = msGMcoarse.vLabels[pv->pvCoarse->idx];
		
		// if we are a right vertex pull value from prolongation table
		else {
			msGMfine.vLabels[pv->idx] = pv->vlblPlongTab[msGMcoarse.vLabels[pv->pvCoarse->idx]];
		}
	}

	// we touched the finer scale's labels
	msGMfine.bLabels = true;
}

//Make Compatible Relaxations graph
template<class GM, class ALG>
void MSMRF<GM, ALG>::makeCRgraph(MSGM& msGMcoarse, MSGM& msGM, MSGM& msGMcr, vector<IndexType>& vidxMapCR2fine)
{

	int numVars = msGM.vLabels.size();
	vector<int> vidxMapFine2CR(numVars, numVars + 1);
	vector<Vertex*> vvrtxFID(numVars);
	IndexType cntr = 0;
	IndexType shp[2] = { msGMcoarse.nLabels, msGMcoarse.nLabels };
	IndexType pidx[2];

	// interpolate labels from coarse grid
	for (int i = 0; i < msGM.vVertices.size(); i++) {
		Vertex*& pv = msGM.vVertices[i];

		// if we are the left vertex in the coarse graph
		if (pv->isLeft()) {
			msGM.vLabels[pv->idx] = msGMcoarse.vLabels[pv->pvCoarse->idx];
		}

		// vertex goes to CR graph
		// map from fine -> CR
		// map from CR -> fine
		else {
			// add variables to CR graph

			msGMcr.vLabels.emplace_back(pv->vlblPlongTab[msGMcoarse.vLabels[pv->pvCoarse->idx]]);	//initialize labels
			vidxMapCR2fine.emplace_back(pv->idx);
			vidxMapFine2CR[pv->idx] = cntr;
			vvrtxFID[cntr] = msGMcr.addVertex(pv,true);
			vvrtxFID[cntr]->setUnary(pv->fn());
			//Vertex* myVertex = msGMcr.addVertex(pv,true);
			cntr++;
		}
	}


	// loop over all edges and load them onto CR graph
	for (int i = 0; i < msGM.vEdges.size(); i++) {
		Edge& edge = *msGM.vEdges[i];

		if ((vidxMapFine2CR[edge.v1.idx] < numVars) && (vidxMapFine2CR[edge.v2.idx] < numVars)) {
			// both ends of the edge are inside the CR graph
			
			// add edge to CR graph
			pidx[0] = vidxMapFine2CR[edge.v1.idx];
			pidx[1] = vidxMapFine2CR[edge.v2.idx];

			msGMcr.addEdge(msGMcr.addFunction(edge.fn()),pidx, pidx + 2);

		}
		else if (vidxMapFine2CR[edge.v1.idx] < numVars) {
			// only left vertex is inside the CR graph

			// load the pairwise onto left vertex' unary
			ExFunction fnUnary(shp, shp + 1);
			LabelType lbl[2];
			lbl[1] = msGM.vLabels[edge.v2.idx];
			for (lbl[0] = 0; lbl[0] < msGM.nLabels; lbl[0]++){
				fnUnary(lbl[0]) = edge(lbl);
			}
			vvrtxFID[vidxMapFine2CR[edge.v1.idx]]->setUnary(fnUnary);

		}
		else if (vidxMapFine2CR[edge.v2.idx] < numVars) {
			// only right vertex is inside to CR graph

			// load the pairwise onto right vertex' unary
			ExFunction fnUnary(shp, shp + 1);
			LabelType lbl[2];
			lbl[0] = msGM.vLabels[edge.v1.idx];
			for (lbl[1] = 0; lbl[1] < msGM.nLabels; lbl[1]++){
				fnUnary(lbl[1]) = edge(lbl);
			}
			vvrtxFID[vidxMapFine2CR[edge.v2.idx]]->setUnary(fnUnary);

		}
		else {
			int n = 1;
		}
	}
}

//Prolong from compatible relaxtions
template<class GM, class ALG>
void MSMRF<GM, ALG>::prolongCR(MSGM& msGMcr, MSGM& msGMfine, vector<IndexType>& vidxMapCR2fine)
{

	for (IndexType i = 0; i < msGMcr.vVertices.size(); i++) {
		Vertex*& pv = msGMcr.vVertices[i];

		msGMfine.vLabels[vidxMapCR2fine[i]] = msGMcr.vLabels[i];

	}

	// we touched the finer scale's labels
	msGMfine.bLabels = true;
}

}