#ifndef POISSON_HPP
#define POISSON_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include "Info.hpp"

#include "../SparseLib++/coord_double.h"
#include "../SparseLib++/comprow_double.h"
#include "../SparseLib++/diagpre_double.h"
#include "../Iml/cg.h"

#include "../Integration/Integration.h"
#include "../LinearAlgebra/LinearAlgebra.h"

#include <list>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#undef MAX
#define MAX(A, B) ((A > B) ? A : B)

class Poisson {
	private:
		// Kernel
		struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
		
		// Triangolazione
		typedef CGAL::Triangulation_vertex_base_with_info_2<PointInfo, K> Vb;
		typedef CGAL::Delaunay_mesh_face_base_2<K, 
			CGAL::Constrained_triangulation_face_base_2<K, 
				CGAL::Triangulation_face_base_with_info_2<FaceInfo , K>
			>
		> Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
		typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
		
		typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
		
		typedef double (*func)(double const &, double const &);
	public:
		typedef CDT::Point Point;
		typedef std::list<Point> PointList;

	public:
		Poisson(PointList const &, func, double const & criteria = 0.);
		~Poisson();
		unsigned int refine(double const &);
		void saveToSVG(char const *, bool const);
	private:
		// LinearAlgebra
		typedef Coord_Mat_double Matrix;
		typedef MV_Vector_double Vector;
		
		// Integrator
		typedef Integrator < Method::NewtonCotes<Geometry::Triangle, 2> >  Integrator2;
		
		CDT cdt_;
		Integrator2 NC;
		
		unsigned int n_point;
		
		func forzante;
	private:
		// Base P1
		static double P1_phi0 (double const & x, double const & y) { return 1 - x -y; }
		static double P1_phi1 (double const & x, double const & y) { return x; }
		static double P1_phi2 (double const & x, double const & y) { return y; }
		
		// Base Clement
		static double C_phi0 (double const & x , double const & y) { return 1.; }
		static double C_phi1 (double const & x , double const & y) { return x; }
		static double C_phi2 (double const & x , double const & y) { return y; }
		
		// Generazione griglia
		void insertNodes(PointList const &, double const &);
		void enumNodes();
		void setBConditions();
		
		// Costruzione matrice
		void makeMatrixTermineNoto(Matrix &, Vector &);
		void printPattern(char const *, Matrix const &);
		
		// Risoluzione del sistema lineare
		void solveSystem(Matrix const &, Vector &, Vector const &);
		
		// Salvataggio soluzione
		void saveSolution(Vector const &);
		
		// Calcolo dei residui
		void calcRes();
		
};

/*
 * Costruttore di default
 */
Poisson::Poisson(PointList const &boundary, func f, double const & criteria) : NC() {
	
	forzante = f;
	
	// Costruzione della griglia
	insertNodes(boundary, criteria);
	setBConditions();
	enumNodes();
	
	// Costruzione matrice termine noto
	Matrix A;
	Vector F;	
	makeMatrixTermineNoto(A, F);
	
	// Soluzione del sistema lineare
	Vector x(n_point, 0.);
	solveSystem(A, x, F);
	
	// Salvataggio della soluzione nella griglia
	saveSolution(x);
}

/*
 * Distruttore
 */
Poisson::~Poisson() {
}

/*
 * Salva in un file SVG
 */
void Poisson::saveToSVG(char const *filename, bool const triangulation=true) {
	int H = 500.;
	int W = 500.;
	std::ofstream f;
	f.open(filename);
	f  << "<?xml version=\"1.0\" standalone=\"yes\"?>\n" \
		<< "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n" \
		<< "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
	// Calcolo altezza e larghezza
	double xmin, xmax, ymin, ymax, umin, umax, rmin, rmax;
	double x, y, u, res;
	int color, r, g ,b;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	xmin = xmax = itV->point().x();
	ymin = ymax = itV->point().y();
	umin = umax = itV->info().value();
	++itV;
	while(itV != cdt_.finite_vertices_end()) {
		x = itV->point().x();
		y = itV->point().y();
		u = itV->info().value();
		if (x < xmin) xmin = x;
		if (x > xmax) xmax = x;
		if (y < ymin) ymin = y;
		if (y > ymax) ymax = y;
		if (u > umax) umax = u;
		if (u < umin) umin = u;
		++itV;
	}
	W = H * (xmax -xmin) / (ymax - ymin);
	// Inizio SVG
	f << "<svg width=\""<<W<<"px\" height=\""<<H<<"px\" " \
	  << "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	// Triangolazione
	CDT::Finite_faces_iterator itF = cdt_.finite_faces_begin();
	rmin = itF->info().res();
	rmax = itF->info().res();
	++itF;
	while(itF != cdt_.finite_faces_end()) {
		res = itF->info().res();
		if (res > rmax) rmax = res;
		if (res < rmin) rmin = res;
		++itF;
	}
	itF = cdt_.finite_faces_begin();
	while(itF != cdt_.finite_faces_end()) {
		// Calcolo dei punti
		double x1 = itF->vertex(0)->point().x(); double y1 = itF->vertex(0)->point().y();
		double x2 = itF->vertex(1)->point().x(); double y2 = itF->vertex(1)->point().y();
		double x3 = itF->vertex(2)->point().x(); double y3 = itF->vertex(2)->point().y();
		x1 = (x1 - xmin) / (xmax - xmin); y1 = 1 - (y1 - ymin) / (ymax - ymin);
		x2 = (x2 - xmin) / (xmax - xmin); y2 = 1 - (y2 - ymin) / (ymax - ymin);
		x3 = (x3 - xmin) / (xmax - xmin); y3 = 1 - (y3 - ymin) / (ymax - ymin);
		x1 = x1 * W; y1 = y1 * H;
		x2 = x2 * W; y2 = y2 * H;
		x3 = x3 * W; y3 = y3 * H;
		
		u = 0.;
		u += itF->vertex(0)->info().value();
		u += itF->vertex(1)->info().value();
		u += itF->vertex(2)->info().value();
		u /= 3;
		res = itF->info().res();
		
		// Calcolo dei colori
		color = 255 * (u - umin ) / (umax - umin);
		// color = 255 * (res - rmin ) / (rmax - rmin);
		r = color;
		g = color;
		b = 0;
		
		// Stampa a file
		if (triangulation)
			f << "<g id=\"triangulation\" fill=\"rgb("<<r<<", "<<g<<", "<<b<<")\" stroke=\"white\" width=\"1\">\n";
		else
			f << "<g id=\"triangulation\" fill=\"rgb("<<r<<", "<<g<<", "<<b<<")\" stroke=\"rgb("<<r<<", "<<g<<", "<<b<<")\" width=\"1\">\n";
		f << "\t<polygon points=\""<<static_cast<int>(x1)<<", "<<static_cast<int>(y1)<<" "\
			<< static_cast<int>(x2)<<", "<<static_cast<int>(y2)<<" "\
			<< static_cast<int>(x3)<<", "<<static_cast<int>(y3)<<"\" />\n";
		f << "</g>\n";
		++itF;
	}
	f  << "</svg>\n";
}

/*
 * Pattern della matrice
 */
void Poisson::printPattern(char const * name, Matrix const & A) {
	std::fstream pbm;
	pbm.open(name, std::ios_base::out);
	
	unsigned int cols = A.dim(0);
	unsigned int rows = A.dim(1);
	
	pbm<<"P1"<<std::endl<<cols<<" "<<rows<<std::endl;
	
	for (unsigned int i = 0 ; i < rows; ++i) {
		for (unsigned int j = 0 ; j < cols; ++j) {
			if (A(i, j) == 0)
				pbm<<"1 ";
			else
				pbm<<"0 ";
			pbm<<std::endl;
		}
	}
}

/*
 * Inserimento dei nodi
 */
void Poisson::insertNodes(PointList const &boundary_, double const & criteria) {
	
	// Criterio mesher
	double minimo = 0.;
	double d;
	double x0, y0, x1, y1;
	
	// Vertex handle list
	typedef std::list<CDT::Vertex_handle> VertexList;
	VertexList Pb_;
	
	// Inserimento punti
	PointList::const_iterator itP;
	
	for (itP = boundary_.begin(); itP != boundary_.end(); ++itP)
		Pb_.push_back(cdt_.insert(*itP));
	
	// Inserimento lati constrained
	VertexList::const_iterator itV1 = Pb_.begin();
	VertexList::const_iterator itV2 = Pb_.begin();
	++itV2;
	
	x0 = (*itV1)->point().x();
	y0 = (*itV1)->point().y();
	x1 = (*itV2)->point().x();
	y1 = (*itV2)->point().y();
	
	minimo = std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	
	while(itV2 != Pb_.end()) {
		cdt_.insert_constraint(*itV1, *itV2);
		
		// lunghezza del lato
		x0 = (*itV1)->point().x();
		y0 = (*itV1)->point().y();
		x1 = (*itV2)->point().x();
		y1 = (*itV2)->point().y();
	
		d = std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
		if (d < minimo)
			minimo = d;
		
		// avanzamento
		++itV1;
		++itV2;
	}
	
	// Ultimo lato
	cdt_.insert_constraint(Pb_.front(), Pb_.back());
	
	x0 = Pb_.front()->point().x();
	y0 = Pb_.front()->point().y();
	x1 = Pb_.back()->point().x();
	y1 = Pb_.back()->point().y();

	d = std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	if (d < minimo)
		minimo = d;
	
	if (criteria != 0.)
		minimo = criteria;
	
	// Impostazione dei nodi delal lista al bordo (da non eliminare)
	CDT::Vertex_circulator cV = cdt_.incident_vertices (cdt_.infinite_vertex());
	CDT::Vertex_circulator begin = cV;
	do {
		cV->info().toBoundary();
		++cV;
	} while(cV != begin);
	
	// Raffinamento
	CGAL::refine_Delaunay_mesh_2(cdt_, Criteria(0.125, minimo));
}

/*
 * Inserisce le condizioni al contorno
 */
void Poisson::setBConditions() {
	CDT::Vertex_circulator cV = cdt_.incident_vertices (cdt_.infinite_vertex());
	CDT::Vertex_circulator begin = cV;
	do {
		cV->info().toDirichlet();
		cV->info().value() = 0.;
		++cV;
	} while(cV != begin);
}

/*
 * Associa ad ogni nodo un indice
 */
void Poisson::enumNodes() {
	unsigned int i = 0u;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	while(itV != cdt_.finite_vertices_end()) {
		if (!itV->info().isDirichlet()) {
			itV->info().index() = i;
			++i;
		}
		++itV;
	}
	n_point = i;
}

/*
 * Assemblaggio della matrice di stiffness e del termine noto
 */
void Poisson::makeMatrixTermineNoto(Matrix &A , Vector &F) {
	// Inserimento valori
	Vector Tmp;
	F.newsize(n_point);
	F = 0.;
	// TODO migliorare l'utilizzo delle matrici sparse
	Tmp.newsize(n_point * n_point);
	Tmp = 0.;
	CDT::Finite_faces_iterator itF = cdt_.finite_faces_begin();
	while(itF != cdt_.finite_faces_end()) {
		
		typedef LinearAlgebra::Vector<double,2> vettore;
		
		vettore p0, p1, p2;
		
		// Coordinate vertici del triangolo
		p0[0] = itF->vertex(0)->point().x(); p0[1] = itF->vertex(0)->point().y();
		p1[0] = itF->vertex(1)->point().x(); p1[1] = itF->vertex(1)->point().y();
		p2[0] = itF->vertex(2)->point().x(); p2[1] = itF->vertex(2)->point().y();
		
		// Geometria dell'integratore
		Integrator2::Geometry g(*itF);
		
		// Area del triangolo
		double detJ = g.det();
		
		// Questo serve per evitare un bug nella libreria CGAL
		detJ = std::abs(detJ);
		if (detJ <= 1.0e-4) {
			++itF;
			continue;
		}
		
		// Determinante dello Jacobiano
		double invdetJ = 1./detJ;
		
		// Indici vertici del triangolo
		int i0 = itF->vertex(0)->info().index();
		int i1 = itF->vertex(1)->info().index();
		int i2 = itF->vertex(2)->info().index();
		
		
	
		vettore  gradphi0(0.), gradphi1(0.), gradphi2(0.);
	    
	
	    // Vettori ausiliari
		gradphi0[0] = p0[1] - p2[1] + p0[2] - p0[0];
		gradphi0[1] = p1[1] - p0[1] + p0[0] - p1[0];
		gradphi1[0] = p2[1] - p0[1];
		gradphi1[1] = p0[1] - p1[1];
		gradphi2[0] = p0[0] - p2[0];
		gradphi2[1] = p1[0] - p0[0];
	
		
		// Condizioni al bordo
		bool b0 = !itF->vertex(0)->info().isDirichlet();
		bool b1 = !itF->vertex(1)->info().isDirichlet();
		bool b2 = !itF->vertex(2)->info().isDirichlet();
		
		// Assemblaggio matrice di stiffness
		
		// Termine di diffusione
		if (b0 && b0) Tmp(n_point*i0+i0) += 0.5*invdetJ*dot(gradphi0,gradphi0);
		if (b0 && b1) Tmp(n_point*i0+i1) += 0.5*invdetJ*dot(gradphi0,gradphi1);
		if (b0 && b2) Tmp(n_point*i0+i2) += 0.5*invdetJ*dot(gradphi0,gradphi2);
		if (b1 && b0) Tmp(n_point*i1+i0) += 0.5*invdetJ*dot(gradphi1,gradphi0);
		if (b1 && b1) Tmp(n_point*i1+i1) += 0.5*invdetJ*dot(gradphi1,gradphi1);
		if (b1 && b2) Tmp(n_point*i1+i2) += 0.5*invdetJ*dot(gradphi1,gradphi2);
		if (b2 && b0) Tmp(n_point*i2+i0) += 0.5*invdetJ*dot(gradphi2,gradphi0);
		if (b2 && b1) Tmp(n_point*i2+i1) += 0.5*invdetJ*dot(gradphi2,gradphi1);
		if (b2 && b2) Tmp(n_point*i2+i2) += 0.5*invdetJ*dot(gradphi2,gradphi2);
		
		// Definizione del dominio dell'integrazione
		NC.domain(g);
		
		// Costruzione termine noto
		if (b0) F[i0] += NC.integrateMul<Integrator2::Transform, Integrator2::NoTransform>(forzante, Poisson::P1_phi0);
		if (b1) F[i1] += NC.integrateMul<Integrator2::Transform, Integrator2::NoTransform>(forzante, Poisson::P1_phi1);
		if (b2) F[i2] += NC.integrateMul<Integrator2::Transform, Integrator2::NoTransform>(forzante, Poisson::P1_phi2);
		
		// Avanzamento dell'iteratore
		++itF;
	}
	unsigned int nz = 0;
	for(int i = 0; i < n_point*n_point; ++i) {
		if (Tmp(i) != 0.)
			++nz;
	}
	A.newsize(n_point, n_point, nz);
	// Inserimento nella matrice
	int j = 0;
	for(int i = 0; i < n_point*n_point; ++i) {
		if(Tmp(i) != 0.) {
			// Valore
			A.val(j) = Tmp(i);
			// Indice colonna
			A.col_ind(j) = i % n_point;
			// Indice riga
			A.row_ind(j) = i / n_point;
			j++;
		}
	}
}

/*
 * Soluzione del sistema
 */

void Poisson::solveSystem(Matrix const & A, Vector & x, Vector const & b) {
	double tol = 1.e-6;
	int maxit = n_point;
	int result;
	CompRow_Mat_double Atmp(A);
	
	DiagPreconditioner_double D(Atmp);
	
	result = CG(Atmp, x, b, D, maxit, tol);
}

/*
 * Salvataggio della soluzione
 */
void Poisson::saveSolution(Vector const & x) {
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	while(itV != cdt_.finite_vertices_end()) {
		if (!itV->info().isDirichlet()) {
			itV->info().value() = x[itV->info().index()];
		}
		++itV;
	}
}

/*
 * Calcolo dei residui sulla faccia
 */
void Poisson::calcRes() {
	
	// Integratore
	Integrator2 integr;
	
	// Iteratori
	CDT::Finite_faces_iterator itF;
	CDT::Finite_vertices_iterator itV;
	
	// Norma H1
	double normH1 = 0.;
	
	/* Passo 1 : calcolo dei gradienti in ogni faccia e hk */
	itF = cdt_.finite_faces_begin();
	while (itF != cdt_.finite_faces_end()) {
		
		typedef LinearAlgebra::Vector<double,2> vettore;
		
		typedef LinearAlgebra::Matrix<double,2,2> matrice;
		
		vettore p0(0.), p1(0.), p2(0.);
		
		// Coordinate vertici del triangolo
	    p0[0] = itF->vertex(0)->point().x(); p0[1] = itF->vertex(0)->point().y();
		p1[0] = itF->vertex(1)->point().x(); p1[1] = itF->vertex(1)->point().y();
		p2[0] = itF->vertex(2)->point().x(); p2[1] = itF->vertex(2)->point().y();
		
		
		// Distanze
		double d01 = norm(p0-p1);
		double d02 = norm(p0-p2);
		double d12 = norm(p1-p2);
		
		// Valori nodali
		double u0 = itF->vertex(0)->info().value();
		double u1 = itF->vertex(1)->info().value();
		double u2 = itF->vertex(2)->info().value();
		
		// Coefficienti del gradiente
		matrice M;
		
		double a01 = (u1 - u0) / d01 / d01;
		double a02 = (u2 - u0) / d02 / d02;
		
		// Gradiente
		itF->info().gradx() = a01 * (x1 - x0) + a02 * (x2 - x0);
		itF->info().grady() = a01 * (y1 - y0) + a02 * (y2 - y0);
		
		// max lati
		itF->info().h() = MAX(MAX(d01, d02), d12);
		
		++itF;
	}
	/* Passo 2 : calcolo dei Clement */
	itV = cdt_.finite_vertices_begin();
	while (itV != cdt_.finite_vertices_end()) {
		// Coordinate del punto
		double x = itV->point().x();
		double y = itV->point().y();
		// Matrice di stiffness
		LinearAlgebra::Matrix<double,3,3> M(0.);
		// Vettore termine noto
		LinearAlgebra::Vector<double,3> bx(0.) , by(0.);
		// Vettore soluzione
		LinearAlgebra::Vector<double,3> ax(0.) , ay(0.);
		// Temporanei
		double t;
		
		// Ciclo sulle facce incidenti al vertice
		CDT::Face_circulator cF, end;
		cF = cdt_.incident_faces(itV);
		end = cF;
		do {
			// Matrice di proiezione locale
			if (!cdt_.is_infinite(cF)) {
				Integrator2::Geometry g(*cF);
				integr.domain(g);
				
				M(0,0) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi0, Poisson::C_phi0);
				M(0,1) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi0, Poisson::C_phi1);
				M(0,2) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi0, Poisson::C_phi2);
				M(1,0) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi1, Poisson::C_phi0);
				M(1,1) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi1, Poisson::C_phi1);
				M(1,2) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi1, Poisson::C_phi2);
				M(2,0) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi2, Poisson::C_phi0);
				M(2,1) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi2, Poisson::C_phi1);
				M(2,2) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(Poisson::C_phi2, Poisson::C_phi2);
				
				bx(0) += cF->info().gradx() * integr.integrate<Integrator2::Transform>(Poisson::C_phi0);
				bx(1) += cF->info().gradx() * integr.integrate<Integrator2::Transform>(Poisson::C_phi1);
				bx(2) += cF->info().gradx() * integr.integrate<Integrator2::Transform>(Poisson::C_phi2);
				
				by(0) += cF->info().grady() * integr.integrate<Integrator2::Transform>(Poisson::C_phi0);
				by(1) += cF->info().grady() * integr.integrate<Integrator2::Transform>(Poisson::C_phi1);
				by(2) += cF->info().grady() * integr.integrate<Integrator2::Transform>(Poisson::C_phi2);
			}
			++cF;
		} while(cF != end);
		
		// Soluzione del sistema
		LinearAlgebra::Solver::LU::solve( M , ax , bx );
		LinearAlgebra::Solver::LU::solve( M , ay , by );
		
		// Gradienti
		itV->info().gradCx() = ax(0) + ax(1) * x + ax(2) * y;
		itV->info().gradCy() = ay(0) + ay(1) * x + ay(2) * y;
		
		++itV;
	}
	/* Passo 3: Norma H1 della soluzione uh */
	itF = cdt_.finite_faces_begin();
	while (itF != cdt_.finite_faces_end()) {
		// Dominio di integrazione
		Integrator2::Geometry g(*itF);
		integr.domain(g);
		
		// Funzione
		LinearAlgebra::Matrix<double, 3, 3> M;
		M(0,0) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi0, Poisson::P1_phi0);
		M(0,1) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi0, Poisson::P1_phi1);
		M(0,2) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi0, Poisson::P1_phi2);
		M(1,0) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi1, Poisson::P1_phi0);
		M(1,1) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi1, Poisson::P1_phi1);
		M(1,2) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi1, Poisson::P1_phi2);
		M(2,0) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi2, Poisson::P1_phi0);
		M(2,1) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi2, Poisson::P1_phi1);
		M(2,2) = integr.integrateMul<Integrator2::NoTransform, Integrator2::NoTransform>(Poisson::P1_phi2, Poisson::P1_phi2);
		
		double u0 = itF->vertex(0)->info().value();
		double u1 = itF->vertex(1)->info().value();
		double u2 = itF->vertex(2)->info().value();
		
		normH1 += (
			u0*u0*M(0,0) + u0*u1*M(0,1) + u0*u2*M(0,2) +
			u1*u0*M(1,0) + u1*u1*M(1,1) + u1*u2*M(1,2) +
			u2*u0*M(2,0) + u2*u1*M(2,1) + u2*u2*M(2,2)
		);
		
		// Gradiente
		normH1 += (itF->info().gradx() * itF->info().gradx()) + (itF->info().grady() * itF->info().grady());
		
		++itF;
	}
	normH1 = std::sqrt(normH1);
	/* Passo 4+5 : calcolo laplaciano + residuo  */
	itF = cdt_.finite_faces_begin();
	while (itF != cdt_.finite_faces_end()) {
		
		// Coordinate vertici del triangolo
		double x0 = itF->vertex(0)->point().x(); double y0 = itF->vertex(0)->point().y();
		double x1 = itF->vertex(1)->point().x(); double y1 = itF->vertex(1)->point().y();
		double x2 = itF->vertex(2)->point().x(); double y2 = itF->vertex(2)->point().y();
		
		// Distanze
		double d01 = std::sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
		double d02 = std::sqrt( (x0-x2)*(x0-x2) + (y0-y2)*(y0-y2) );
		double d12 = std::sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
		
		// Valori nodali
		double gx0 = itF->vertex(0)->info().gradCx();
		double gx1 = itF->vertex(1)->info().gradCx();
		double gx2 = itF->vertex(2)->info().gradCx();
		double gy0 = itF->vertex(0)->info().gradCy();
		double gy1 = itF->vertex(1)->info().gradCy();
		double gy2 = itF->vertex(2)->info().gradCy();
		
		// Coefficienti del gradiente
		double ax01 = (gx1 - gx0) / d01;
		double ax02 = (gx2 - gx0) / d02;
		double ay01 = (gy1 - gy0) / d01;
		double ay02 = (gy2 - gy0) / d02;
		
		// Laplaciano
		double lapo = ax01 * (x1 - x0) / d01 + ax02 * (x2 - x0) / d02 +
			ay01 * (y1 - y0) / d01 + ay02 * (y2 - y0) / d02;
		
		// Calcolo del residuo
		itF->info().res() = 0.;
		Integrator2::Geometry g(*itF);
		integr.domain(g);
		
		// Residuo sulla faccia
		itF->info().res() += itF->info().h() * std::sqrt(
			integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(forzante , forzante) + 
			2 * lapo * integr.integrate<Integrator2::Transform>(forzante) +
			lapo * lapo * std::abs(g.det()) * 0.5 
		);
		
		// Residuo al bordo
		double t = 0.;
		
		if (!cdt_.is_infinite(itF->neighbor(0))) {
			t += pow(
				(itF->info().gradx() - itF->neighbor(0)->info().gradx()) * (y1 - y2) +
				(itF->info().grady() - itF->neighbor(0)->info().grady()) * (x2 - x1)
			, 2) / d12;
		}
		
		if (!cdt_.is_infinite(itF->neighbor(1))) {
			t += pow(
				(itF->info().gradx() - itF->neighbor(1)->info().gradx()) * (y2 - y0) +
				(itF->info().grady() - itF->neighbor(1)->info().grady()) * (x0 - x2)
			, 2) / d02;
		}
		
		if (!cdt_.is_infinite(itF->neighbor(2))) {
			t += pow(
				(itF->info().gradx() - itF->neighbor(2)->info().gradx()) * (y0 - y1) +
				(itF->info().grady() - itF->neighbor(2)->info().grady()) * (x1 - x0)
			, 2) / d01;
		}
		
		itF->info().res() += std::sqrt(itF->info().h() * t) * 0.5;
		
		itF->info().res() /= normH1;
		
		++itF;
	}
}

/*
 * Raffinamento della mesh
 */
unsigned int Poisson::refine(double const & eps) {
	
	double max_res = (3. * eps ) / (2. * std::sqrt(cdt_.number_of_faces ()));
	double min_res = max_res / 3.;
	
	PointList da_inserire;
	std::list<CDT::Vertex_handle> da_eliminare;
	
	CDT::Finite_faces_iterator itF;
	CDT::Finite_vertices_iterator itV;
	
	// Calcolo residuo
	calcRes();
	
	// Inserimento
	itF = cdt_.finite_faces_begin();
	while (itF != cdt_.finite_faces_end()) {
		Integrator2::Geometry g(*itF);
		if ((itF->info().res() > max_res) && (std::abs(g.det()) > 9.e-4)) {
			double x0 = itF->vertex(0)->point().x(); double y0 = itF->vertex(0)->point().y();
			double x1 = itF->vertex(1)->point().x(); double y1 = itF->vertex(1)->point().y();
			double x2 = itF->vertex(2)->point().x(); double y2 = itF->vertex(2)->point().y();
			da_inserire.push_back(Point(
				(x0+x1+x2)/3.,
				(y0+y1+y2)/3.
			));
		}
		++itF;
	}
	
	// Rimozione
	itV =cdt_.finite_vertices_begin();
	while (itV != cdt_.finite_vertices_end()) {
		
		// Non eliminabili
		if (itV->info().isBoundary()) {
			++itV;
			continue;
		}
		
		// Residuo medio
		double t = 0.;
		unsigned int n = 0;
		CDT::Face_circulator cF = cdt_.incident_faces(itV);
		CDT::Face_circulator end = cF;
		do {
			if (!cdt_.is_infinite(cF)) {
				t += cF->info().res();
				++n;
			}
			++cF;
		} while(cF != end);
		t /= n;
		
		// Inserimento
		if (t < min_res)
			da_eliminare.push_back(itV);
		
		++itV;
	}
	
	if (da_eliminare.size()) {
		std::list<CDT::Vertex_handle>::iterator it = da_eliminare.begin();
		while(it != da_eliminare.end()) {
			cdt_.remove(*it);
			++it;
		} 
	}
	if (da_inserire.size())
		cdt_.insert(da_inserire.begin(), da_inserire.end());
	
	// Costruzione della griglia
	setBConditions();
	enumNodes();
	
	// Costruzione matrice termine noto
	Matrix A;
	Vector F;	
	makeMatrixTermineNoto(A, F);
	
	// Soluzione del sistema lineare
	Vector x(n_point, 0.);
	solveSystem(A, x, F);
	
	// Salvataggio della soluzione nella griglia
	saveSolution(x);
	
	return da_inserire.size() + da_eliminare.size();
}

#endif // POISSON_HPP
