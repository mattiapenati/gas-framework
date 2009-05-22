#ifndef INFO_HPP
#define INFO_HPP

/* informazioni sui punti */
class PointInfo {
	public:
		PointInfo();
		~PointInfo();
		PointInfo &operator=(PointInfo const &);
		unsigned int & index();
		double & value();
		
		void toDirichlet();
		void toNeumann();
		
		void toBoundary();
		
		bool isDirichlet();
		
		bool isBoundary();
		
	private:
		enum BConditionType { NONE, DIRICHLET, NEUMANN };
		
		unsigned int index_;
		double value_;
		BConditionType bcond_;
		bool bound_list_;
};

/* informazioni sulle facce */
class FaceInfo {
	public:
		FaceInfo();
		~FaceInfo();
		FaceInfo & operator= (FaceInfo const &);
		
		double & res();
	
	private:
		/* residuo (stimatore H1) */
		double res_;
};

PointInfo::PointInfo(): index_(0u), value_(), bcond_(NONE), bound_list_(false) {
}

PointInfo::~PointInfo() {
}

PointInfo & PointInfo::operator=(PointInfo const &i) {
	index_ = i.index_;
	value_ = i.value_;
	bcond_ = i.bcond_;
	bound_list_ = i.bound_list_;
	
	return *this;
}

unsigned int &PointInfo::index() {
	return index_;
}

double &PointInfo::value() {
	return value_;
}

void PointInfo::toDirichlet() {
	bcond_ = DIRICHLET;
}

void PointInfo::toNeumann() {
	bcond_ = NEUMANN;
}

void PointInfo::toBoundary() {
	bound_list_ = true;
}

bool PointInfo::isDirichlet() {
	return (bcond_ == DIRICHLET);
}

bool PointInfo::isBoundary() {
	return bound_list_;
}

FaceInfo::FaceInfo(): res_() {
}

FaceInfo::~FaceInfo() {
}

FaceInfo & FaceInfo::operator= (FaceInfo const & f) {
	res_ = f.res_;
}

double & FaceInfo::res() {
	return res_;
}

#endif // INFO_HPP