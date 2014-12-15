#pragma once
#include "MultSclVertex.hxx"
#include "MultSclFactor.hxx"

using namespace std;

namespace opengm {

template<class GM>
class MultSclGM :
	public GM
{
public:
	// get the OpenGM typedefs
	typedef GM GraphicalModelType;
	OPENGM_GM_TYPE_TYPEDEFS;
	typedef ExplicitFunction<double,LabelType,LabelType> ExFunction;

	// <summary> Default constructor.</summary>
	// <remarks> Constructs an empty multiscale grpahical layer.</remarks>
	MultSclGM();

	// <summary> Constructs from an existing graphical model.</summary>
	// <remarks> Used to construct a multiscale graphical layer from an already existing graphical
	// model.</remarks>
	// <param name="gm"> The existing single scale graphical model.</param> 
	MultSclGM(GM& gm);

	MultSclGM(IndexType nLabels);

	// <summary> Destructor.</summary>
	~MultSclGM(void);

	vector<LabelType> vLabels;		///< The current state's labels
	bool bLabels;					///< flagged when labels are firstly assigned
	IndexType nLabels;				///< Total number of labels in the model

	typedef MultSclFactor<MultSclGM<GM>> Edge;		///< Multiscale factor representing an edge
	vector<Edge*> vEdges;		///< Vector of edges (dual-var factors)
	typedef MultSclVertex<MultSclGM<GM>> Vertex;	///< Multiscale vertex representing a variable
	vector<Vertex*> vVertices;	///< Multiscale vertices representing the variables

	inline void reserveEdges(int nEdges);
	inline void reserveVertices(int nVertices);
	Edge& addEdge(const FunctionIdentifier & fnIdent, LabelType* plblBegin, LabelType* plblEnd);
	Vertex& addVertex(Vertex* pvFine);
	Vertex* addVertex(Vertex* pvFine, bool bPtr);
	

private:
		
	MultSclGM(MultSclGM& msGM);	//deleted

};


template<class GM>
MultSclGM<GM>::MultSclGM() :
	GM()
{
}

template<class GM>
MultSclGM<GM>::~MultSclGM()
{
	for_each(vEdges.begin(), vEdges.end(), [](Edge*& pEdge) { delete pEdge; });
	for_each(vVertices.begin(), vVertices.end(), [](Vertex*& pv) { delete pv; });
}

template<class GM>
MultSclGM<GM>::MultSclGM(GM& gm) :
	GM(),	// construct using the base class copy constructor
	nLabels(gm.numberOfLabels(0)),
	vLabels(gm.numberOfVariables()),
	bLabels(false)
{
	// reserve space for our vertices and edges
	vVertices.resize(gm.numberOfVariables());
	vEdges.reserve(gm.numberOfFactors() - gm.numberOfVariables());
	
#ifdef DEBUG
	/*int n = gm.numberOfVariables();
	int ff = gm.numberOfFactors();
	vector<LabelType> vlbl(gm.numberOfVariables());
	for (int v= 0; v < vlbl.size(); v++)
		vlbl[v] = (LabelType)floor((double)rand() / (double)RAND_MAX * 3 + 0.5);
	double d = gm.evaluate(vlbl);*/
#endif

	for (IndexType v = 0; v < gm.numberOfVariables(); v++)
		addVariable(nLabels);

	IndexType shp[2] = {nLabels, nLabels};
	map<FunctionIdentifier, FunctionIdentifier> mapFI;

	// iterate unary factors and create vertices
	for (IndexType f = 0, iFactor; f < gm.numberOfFactors(); f++) {
		const typename GM::FactorType& fac = gm[f];
		FunctionIdentifier fiOld(fac.functionIndex(), fac.functionType());

		if (gm.numberOfVariables(f) == 1) {
#ifdef DEBUG
			int nVarLabels = gm.numberOfLabels(fac.variableIndex(0));
			//assert(nVarLabels == nLabels);
#endif
			ExFunction fn(shp, shp + 1);
			for (LabelType lbl = 0; lbl < nLabels; lbl++)
				fn(lbl) = fac(&lbl);

			iFactor = addFactor(addFunction(fn),
				fac.variableIndicesBegin(), fac.variableIndicesEnd());
			vVertices[fac.variableIndex(0)] = (new Vertex(*this, fac.variableIndex(0), iFactor));
		}
	}
	//iterate variables which don't have a unary and create them with zero cost
	for (IndexType v = 0, iFactor; v < vVertices.size(); v++) {
		if (!vVertices[v]) {
			ExFunction fn(shp, shp + 1);
			for (LabelType lbl = 0; lbl < nLabels; lbl++)
				fn(lbl) = 0;

			iFactor = addFactor(addFunction(fn), &v, &v + 1);
			vVertices[v] = (new Vertex(*this, v, iFactor));
		}
	}
	//iterate binary factors and create edges
	for (IndexType f = 0, iFactor; f < gm.numberOfFactors(); f++) {
		const typename GM::FactorType& fac = gm[f];
		FunctionIdentifier fiOld(fac.functionIndex(), fac.functionType());

		if (gm.numberOfVariables(f) == 2) {
			typename map<FunctionIdentifier, FunctionIdentifier>::iterator 
				itr = mapFI.find(fiOld);
			FunctionIdentifier fiNew;
			if (itr == mapFI.end()) {
				ExFunction fn(shp, shp+2);
				LabelType lbl[2];
				for (lbl[0] = 0; lbl[0] < nLabels; lbl[0]++) {
					for (lbl[1] = 0; lbl[1] < nLabels; lbl[1]++)
						fn(lbl) = fac(lbl);
				}
				
				fiNew = addFunction(fn);
				mapFI[fiOld] = fiNew;
			}
			else
				fiNew = itr->second;

#ifdef DEBUG
			/*VectorView<vector<long unsigned int>, long unsigned int>::const_iterator p = fac.variableIndicesBegin();;
			IndexType vv[10];
			for (int i = 0; p != fac.variableIndicesEnd(); i++) {
				vv[i] = *p;
				p++;
			}*/
#endif 

			iFactor = addFactor(fiNew, fac.variableIndicesBegin(), fac.variableIndicesEnd());
			vEdges.emplace_back(new Edge(*this, iFactor));
		}
	}

#ifdef DEBUG
	/*double dNew = evaluate(vlbl);*/
#endif

}

template<class GM>
MultSclGM<GM>::MultSclGM(IndexType nLabels) :
	GM(), nLabels(nLabels), bLabels(false)
{
}

template<class GM>
inline void MultSclGM<GM>::reserveEdges(int nEdges) 
{
	vEdges.reserve(nEdges);
}

template<class GM>
inline void MultSclGM<GM>::reserveVertices(int nVertices) 
{
	vVertices.reserve(nVertices);
	vLabels.reserve(nVertices);
}

template<class GM>
typename MultSclGM<GM>::Edge& MultSclGM<GM>::addEdge(const FunctionIdentifier& fnIdent, 
								   LabelType* plblBegin, LabelType* plblEnd)
{
	IndexType iFactor = addFactor(fnIdent, plblBegin, plblEnd);
	vEdges.emplace_back(new Edge( *this, iFactor ));

	return *vEdges.back();
}

template<class GM>
typename MultSclGM<GM>::Vertex& MultSclGM<GM>::addVertex(Vertex* pvFine)
{
	// create a coarse vertex as the parent of the supplied fine one
	vVertices.emplace_back(new Vertex(*this, vVertices.size(), pvFine));
	Vertex& v = *vVertices.back();
	
	// add the variable to the model
	addVariable(nLabels);

	return v;
}

template<class GM>
typename MultSclGM<GM>::Vertex* MultSclGM<GM>::addVertex(Vertex* pvFine, bool bPtr)
{
	// create a coarse vertex as the parent of the supplied fine one
	vVertices.emplace_back(new Vertex(*this, vVertices.size(), pvFine));
	Vertex* v = vVertices.back();

	// add the variable to the model
	addVariable(nLabels);

	return v;
}

}
