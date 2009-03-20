#ifndef INFO_HPP
#define INFO_HPP

class PointInfo {
	public:
		PointInfo();
		~PointInfo();
		PointInfo &operator=(PointInfo const &);
		unsigned int &index();
		double &value();
		
		void toDirichlet();
		void toNeumann();
		
		bool isDirichlet();
		
	private:
		enum BConditionType { NONE, DIRICHLET, NEUMANN };
		
		unsigned int index_;
		double value_;
		BConditionType bcond_;
};

PointInfo::PointInfo(): index_(0u), value_(0.), bcond_(NONE) {
}

PointInfo::~PointInfo() {
}

PointInfo &PointInfo::operator=(PointInfo const &i) {
	index_ = i.index_;
	value_ = i.value_;
	bcond_ = i.bcond_;
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

bool PointInfo::isDirichlet() {
	return (bcond_ == DIRICHLET);
}

#endif // INFO_HPP