#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <list>
#include <ostream>

class Info {
	public:
		Info();
		Info(Info &);
		~Info();
		Info &operator=(Info &);
		// Condizioni
		void toDirichlet();
		void toNeumann();
		void toNone();
		// Campi
		double &value();
		unsigned int &index();
	private:
		enum BoundaryCondition {
			NONE = 0,
			DIRICHLET = 1,
			NEUMANN = 2
		};
		
		double value_;
		BoundaryCondition condition_;
		unsigned int index_;
};

class Poisson {
	private:
		// Kernel
		struct K: CGAL::Exact_predicates_inexact_constructions_kernel {};
		
		// Struttura dati
		typedef CGAL::Triangulation_vertex_base_with_info_2<Info,K> Vb; 
		typedef CGAL::Triangulation_face_base_2<K> Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
		
		// Triangolazione
		typedef CGAL::Constrained_Delaunay_triangulation_2<K,Tds> CDT;
		
		// Algebra lineare
		typedef CompRow_Mat_double Matrix;
		typedef MV_Vector_double Vector;
	public:
		typedef CDT::Point Point;
		typedef std::list<Point> PointList;
		
		Poisson(PointList &);
		~Poisson();
		void solve();
		void saveToFile(std::ostream &);
		void saveToFile(std::ostream &, std::ostream &);
		void saveToSvg(std::	ostream &);
	private:
		CDT cdt_;
		PointList boundary_;
		
		// Costruisce la triangolazione
		void makeTriangulation_();
		// Costruisce la matrice di stiffness
		void makeMatrix_(Matrix &);
		// Costruisce il vettore del termine noto
		void makeF_(Vector &);
		// Risolve il sistema lineare
		void solve_(Matrix &, Vector &, Vector &);
		// Copia la soluzione nella tringolazione
		void saveSolution_(Vector &);
};
