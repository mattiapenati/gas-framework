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
		
		double & gradCx();
		double & gradCy();
		
	private:
		enum BConditionType { NONE, DIRICHLET, NEUMANN };
		
		unsigned int index_;
		double value_;
		BConditionType bcond_;
		bool bound_list_;
		
		/* Gradiente della proiezione di Clement */
		double gradCx_;
		double gradCy_;
};

/* informazioni sulle facce */
class FaceInfo {
	public:
		FaceInfo();
		~FaceInfo();
		FaceInfo & operator= (FaceInfo const &);
		
		double & gradx();
		double & grady();
		
		double & h();
		
		double & res();
	
	private:
		/* gradiente della soluzione uh (sui P1) */
		double gradx_;
		double grady_;
		/* ampiezza della faccia */
		double h_;
		/* residuo (stimatore H1) */
		double res_;
};

PointInfo::PointInfo(): index_(0u), value_(), bcond_(NONE), gradCx_(), gradCy_(), bound_list_(false) {
}

PointInfo::~PointInfo() {
}

PointInfo & PointInfo::operator=(PointInfo const &i) {
	index_ = i.index_;
	value_ = i.value_;
	bcond_ = i.bcond_;
	gradCx_ = i.gradCx_;
	gradCy_ = i.gradCy_;
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

double & PointInfo::gradCx() {
	return gradCx_;
}
double & PointInfo::gradCy() {
	return gradCy_;
}

FaceInfo::FaceInfo(): gradx_(), grady_(), h_(), res_() {
}

FaceInfo::~FaceInfo() {
}

FaceInfo & FaceInfo::operator= (FaceInfo const & f) {
	gradx_ = f.gradx_;
	grady_ = f.grady_;
	h_ = f.h_;
	res_ = f.res_;
}

double & FaceInfo::gradx() {
	return gradx_;
}
double & FaceInfo::grady() {
	return grady_;
}
double & FaceInfo::h() {
	return h_;
}
double & FaceInfo::res() {
	return res_;
}

#endif // INFO_HPP