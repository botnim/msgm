// file:	MultSclVertex.hxx
//
// summary:	Definition and implementation of the multi-scale vertex class

#pragma once
#include "MultSclFactor.hxx"

using namespace std;

namespace opengm {

// <summary> MultSclVertex class - A multi-scale vertex.</summary>
// <remarks> This is one of the "helper" classes to the multiscale framework. It represents a
//  variable ("vertex") in the model and includes extra information on it such as the pointer to
//  the interpolating vertex from the finer scale, an array of interpolant variables with their
//  enhanced factors, and an array of current prior statistical info for the final state .</remarks>
// 
// <typeparam name="MSGM"> The type multiscale grpahical model this vertex corresponds to.</typeparam>

template<class MSGM>
class MultSclVertex
{
private:

	// get the OpenGM typedefs
	typedef MSGM GraphicalModelType;
	OPENGM_GM_TYPE_TYPEDEFS;
	typedef typename MSGM::ExFunction ExFunction;
	typedef typename MSGM::Vertex Vertex;


public:

	MultSclVertex() { assert(false); };	//deleted
	MultSclVertex(const MultSclVertex<MSGM>& v) : 
		msGM(v.msGM), iFactor(v.iFactor), fi(v.fi) { assert(false); }; //deleted

	// <summary> Construct a multiscale vertex.</summary>
	// <remarks> Construct a multiscale vertex and binds it to the specified multiscale graphical
	//  model.</remarks>
	// <param name="msGM">    [in,out] The multiscale graphical model this vertex is bound to.</param>
	// <param name="idx">	  The index.</param>
	// <param name="pFactor"> [in,out] If non-null, the factor.</param>
	MultSclVertex(MSGM& msGM, IndexType idx, IndexType iFactor);

	// <summary> Constructor.</summary>
	// <remarks> Stav, 16/01/2014.</remarks>
	// <param name="msGM">   [in,out] The milliseconds gm.</param>
	// <param name="idx">    The index.</param>
	// <param name="pvLeft"> [in,out] If non-null, the pv left.</param>
	MultSclVertex(MSGM& msGM, IndexType idx, Vertex* pvLeft);

	~MultSclVertex(void);
	
private:

	MSGM& msGM;									///< The multiscale graphical model this vertex belongs to
	IndexType iFactor;							///< Unary factor (index) for this vertex
	vector<double> vdPostPred;					///< Posterior predictive value for each label
	double dPostPredSum;						///< Simply holds the sum of the above
	
public:

	FunctionIdentifier fi;						///< Function identifier for unary term
	typename MSGM::Vertex* pvCoarse;			///< Pointer to the "coarse" vertex this vertex belongs to
	typename MSGM::Vertex* pvLeft;				///< Pointer to the "left" vertex
	vector<typename MSGM::Vertex*> vpvRights;	///< Pointers to "right" vertices
	//vector<typename MSGM::Edge*> vpEdges;		///< All the edges ("enhanced") included in the vertex
	
	vector<double>* pvdPostPred;				
	double* pdPostPredSum;						
	IndexType idx;						///< The index

	vector<LabelType> vlblPlongTab;		///< The prolongation table

	ValueType operator()(LabelType lbl);
	typename MSGM::ExFunction& fn();
	void addRight(Vertex* pvRight);
	bool isLeft();
	void setUnary(typename MSGM::ExFunction& fnUnary);
	
};

// construct a fine level vertex
template<class MSGM>
MultSclVertex<MSGM>::MultSclVertex(MSGM& msGM, IndexType idx, IndexType iFactor) :
	msGM(msGM), idx(idx), iFactor(iFactor),
	fi(FunctionIdentifier(msGM[iFactor].functionIndex(), msGM[iFactor].functionType())),
	pvCoarse(NULL),
	pvLeft(NULL),
	vdPostPred(msGM.nLabels, 1),		// initialize post predictive with a uniform 1 for
	dPostPredSum((double)msGM.nLabels),	//	for every label and also the appropriate sum
	pvdPostPred(&vdPostPred),
	pdPostPredSum(&dPostPredSum),
	vlblPlongTab(msGM.nLabels)
{
}

// construct a coarse vertex
template<class MSGM>
MultSclVertex<MSGM>::MultSclVertex(MSGM& msGM, IndexType idx, Vertex* pvLeft) :
	msGM(msGM), idx(idx), pvLeft(pvLeft), pvCoarse(NULL), 
	pvdPostPred(pvLeft->pvdPostPred),	//point the post-pred values to the fine left vertex's
	pdPostPredSum(pvLeft->pdPostPredSum),
	vlblPlongTab(msGM.nLabels)
{
	//mainting linking with left child
	pvLeft->pvCoarse = this;
}

template<class MSGM>
MultSclVertex<MSGM>::~MultSclVertex()
{
	// coarse vertices get destroyed before their "children", we therefor free descendant
	// vertices from the links to the current destroyed one. 
	if (pvLeft) {
		pvLeft->pvCoarse = NULL;
		for_each(vpvRights.begin(), vpvRights.end(), 
			[](Vertex*& pv) { pv->pvCoarse = NULL; });
	}
}

template<class MSGM>
typename MultSclVertex<MSGM>::ValueType MultSclVertex<MSGM>::operator()(LabelType lbl)
{
	return msGM[iFactor](&lbl);
}

template<class MSGM>
void MultSclVertex<MSGM>::addRight(Vertex* pvRight)
{
	//add to rights vector
	vpvRights.push_back(pvRight);

	//maintain linking
	pvRight->pvCoarse = this;
}

template<class MSGM>
bool MultSclVertex<MSGM>::isLeft()
{
	return pvCoarse->pvLeft == this;
}

template<class MSGM>
void MultSclVertex<MSGM>::setUnary(typename MSGM::ExFunction& fnUnary)
{
	fi = msGM.addFunction(fnUnary);
	iFactor = msGM.addFactor(fi, &idx, &idx+1);
}

template<class MSGM>
inline typename MSGM::ExFunction& MultSclVertex<MSGM>::fn()
{
	return msGM.template getFunction<ExFunction>(fi);
}

}