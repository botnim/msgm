#pragma once

using namespace std;

namespace opengm {

template<typename MSGM>
class MultSclFactor
{
public:
	// get the OpenGM typedefs
	typedef MSGM GraphicalModelType;
	OPENGM_GM_TYPE_TYPEDEFS;

	MultSclFactor(void) { assert(false); };	//deleted;
	MultSclFactor(const MultSclFactor<MSGM>& f) : msGM(f.msGM), 
		iFactor(f.iFactor), v1(f.v1), v2(f.v2), fi(f.fi) { assert(false); };	//deleted

	MultSclFactor(MSGM& msGM, IndexType iFactor);
	~MultSclFactor(void);

private:

	typedef typename MSGM::ExFunction ExFunction;
	
	MSGM& msGM;						///< The multiscale graphical model this vertex belongs to
	IndexType iFactor;
	
public:

	FunctionIdentifier fi;			///< Function identifier for unary term
	typename MSGM::Vertex &v1,&v2;	///< the vertices
	size_t dwUser;  ///< Used to assocciate different data with the edge

	ValueType operator()(LabelType* begin, bool bRev);
	ValueType operator()(LabelType* begin);
	ExFunction& fn();
};


template<class MSGM>
MultSclFactor<MSGM>::MultSclFactor(MSGM& msGM, IndexType iFactor) :
	msGM(msGM), iFactor(iFactor),
	fi(FunctionIdentifier(msGM[iFactor].functionIndex(), msGM[iFactor].functionType())),
	/*v1(msGM[iFactor].variableIndex(0) < msGM[iFactor].variableIndex(1) ? 
		*msGM.vVertices[msGM[iFactor].variableIndex(0)] : *msGM.vVertices[msGM[iFactor].variableIndex(1)]),
	v2(msGM[iFactor].variableIndex(0) < msGM[iFactor].variableIndex(1) ? 
		*msGM.vVertices[msGM[iFactor].variableIndex(1)] : *msGM.vVertices[msGM[iFactor].variableIndex(0)])*/
	v1(*msGM.vVertices[msGM[iFactor].variableIndex(0)]),
	v2(*msGM.vVertices[msGM[iFactor].variableIndex(1)])
{
#ifdef DEBUG
	assert(&v1 != 0 && &v2 != 0);
#endif
}

template<class MSGM>
MultSclFactor<MSGM>::~MultSclFactor()
{
}

template<class MSGM>
typename MultSclFactor<MSGM>::ValueType MultSclFactor<MSGM>::operator()(LabelType* begin, bool bRev)
{
	//check if edge is "transposed"
	if (!bRev) {
		return msGM[iFactor](begin);
	}
	else {
		LabelType lbl_[2];
		lbl_[0] = begin[1];
		lbl_[1] = begin[0];
		return msGM[iFactor](lbl_);
	}
}

template<class MSGM>
typename MultSclFactor<MSGM>::ValueType MultSclFactor<MSGM>::operator()(LabelType* begin)
{
	return msGM[iFactor](begin);
}

template<class MSGM>
inline typename MSGM::ExFunction& MultSclFactor<MSGM>::fn()
{
	return msGM.template getFunction<ExFunction>(fi);
}


}
