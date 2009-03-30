#ifndef POISSON_HPP
#define POISSON_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include "Info.hpp"

#include "Gas/SparseLib++/comprow_double.h"

#include "Gas/SparseLib++/mv/mvvi.h"
#include <list>
#include <fstream>
#include <iostream>
#include <string>

class Poisson {
	private:
		// Kernel
		struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
		
		// Triangolazione
		typedef CGAL::Triangulation_vertex_base_with_info_2<PointInfo, K> Vb;
		typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
		typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
		
		typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
	public:
		typedef CDT::Point Point;
		typedef std::list<Point> PointList;
	public:
		Poisson(PointList const &);
		~Poisson();
		void saveToSVG(char const *, bool const);
	private:
		// LinearAlgebra
		typedef CompRow_Mat_double Matrix;
		typedef MV_Vector_double Vector;		
		typedef MV_Vector_int VectorInt;

		CDT cdt_;
	private:
		typedef std::list<CDT::Vertex_handle> VertexList;
		
		// Generazione griglia
		void insertNodes(PointList const &);
		void enumNodes();
		void setBConditions();
		// Costruzione matrice
		void makeMatrix(Matrix &);
		void makeTermineNoto(Vector &);
};

/*
 * Costruttore di default
 */
Poisson::Poisson(PointList const &boundary) {
	// Costruzione della griglia
	insertNodes(boundary);
	setBConditions();
	enumNodes();
	// Costruzione matrice
	Matrix A;
	makeMatrix(A);
	// Costruzione termine noto
	Vector F;
	makeTermineNoto(F);
	// Soluzione del sistema lineare
	// Salvataggio della soluzione nella griglia
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
	double xmin, xmax, ymin, ymax;
	double x, y;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	xmin = xmax = itV->point().x();
	ymin = ymax = itV->point().y();
	++itV;
	while(itV != cdt_.finite_vertices_end()) {
		x = itV->point().x();
		y = itV->point().y();
		if (x < xmin) xmin = x;
		if (x > xmax) xmax = x;
		if (y < ymin) ymin = y;
		if (y > ymax) ymax = y;
		++itV;
	}
	W = H * (xmax -xmin) / (ymax - ymin);
	// Inizio SVG
	f  << "<svg width=\""<<W<<"px\" height=\""<<H<<"px\" " \
		<< "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	// Triangolazione
	if (triangulation) {
		f << "<g id=\"triangulation\" fill=\"none\" stroke=\"black\" widht=\"1\">\n";
		CDT::Finite_faces_iterator itF = cdt_.finite_faces_begin();
		while(itF != cdt_.finite_faces_end()) {
			double x1 = itF->vertex(0)->point().x(); double y1 = itF->vertex(0)->point().y();
			double x2 = itF->vertex(1)->point().x(); double y2 = itF->vertex(1)->point().y();
			double x3 = itF->vertex(2)->point().x(); double y3 = itF->vertex(2)->point().y();
			x1 = (x1 - xmin) / (xmax - xmin); y1 = 1 - (y1 - ymin) / (ymax - ymin);
			x2 = (x2 - xmin) / (xmax - xmin); y2 = 1 - (y2 - ymin) / (ymax - ymin);
			x3 = (x3 - xmin) / (xmax - xmin); y3 = 1 - (y3 - ymin) / (ymax - ymin);
			x1 = x1 * W; y1 = y1 * H;
			x2 = x2 * W; y2 = y2 * H;
			x3 = x3 * W; y3 = y3 * H;
			f  << "<polygon points=\""<<static_cast<int>(x1)<<", "<<static_cast<int>(y1)<<" "\
				<< static_cast<int>(x2)<<", "<<static_cast<int>(y2)<<" "\
				<< static_cast<int>(x3)<<", "<<static_cast<int>(y3)<<"\" />\n";
			++itF;
		}
		f << "</g>\n";
	}
	f  << "</svg>\n";
}

/*
 * Inserimento dei nodi
 */
void Poisson::insertNodes(PointList const &boundary_) {
	// Vertex handle list
	VertexList Pb_;
	
	// Inserimento punti
	PointList::const_iterator itP;
	for (itP = boundary_.begin(); itP != boundary_.end(); ++itP)
		Pb_.push_back(cdt_.insert(*itP));
	
	// Inserimento lati constrained
	VertexList::const_iterator itV1 = Pb_.begin();
	VertexList::const_iterator itV2 = Pb_.begin();
	++itV2;
	while(itV2 != Pb_.end()) {
		cdt_.insert_constraint(*itV1, *itV2);
		++itV1;
		++itV2;
	}
	
	// Ultimo lato
	cdt_.insert_constraint(Pb_.front(), Pb_.back());
	
	// Raffinamento
	CGAL::refine_Delaunay_mesh_2(cdt_, Criteria(0.125, 0.5));
}

/*
 * Inserisce le condizioni al contorno
 */
void Poisson::setBConditions() {
	CDT::Vertex_circulator cV = cdt_.incident_vertices (cdt_.infinite_vertex());
	CDT::Vertex_circulator begin = cV;
	do {
		cV->info().toDirichlet();
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
}

/*
 * Assemblaggio della matrice di stiffness
 */
void Poisson::makeMatrix(Matrix &A) {
	// Calcolo degli elementi non nulli e dimensione del problema
	unsigned int nz = 0u;
	unsigned int n = 0u;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	CDT::Vertex_circulator cV;
	CDT::Vertex_circulator begin;
	while (itV != cdt_.finite_vertices_end()) {
		cV = cdt_.incident_vertices(itV);
		if (!(itV->info().isDirichlet()))
			++n;
		begin = cV;
		do {
			if (!cdt_.is_infinite(cV))
				++nz;
			++cV;
		} while(cV != begin);
		++itV;
	}
	// Inserimento valori
	A.newsize(n, n, nz);
	CDT::Finite_faces_iterator itF = cdt_.finite_faces_begin();
	Vector Tmp;
	Tmp.newsize(n*n);
	Tmp = 0;
	while (itF != cdt_.finite_faces_end()) {
		// Coordinate vertici del triangolo
		double x0 = itF->vertex(0)->point().x();
		double y0 = itF->vertex(0)->point().y();
		double x1 = itF->vertex(1)->point().x();
		double y1 = itF->vertex(1)->point().y();
		double x2 = itF->vertex(2)->point().x();
		double y2 = itF->vertex(2)->point().y();
		// Indici vertici del triangolo
		int i0 = itF->vertex(0)->info().index();
		int i1 = itF->vertex(1)->info().index();
		int i2 = itF->vertex(2)->info().index();
		// Determinante dello Jacobiano
		double invdetJ = 1/((x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0));
		// Vettori ausiliari
		double phi01 = y0 - y2 + x2 - x0;
		double phi02 = y1 - y0 + x0 - x1;
		double phi11 = y2 - y0;
		double phi12 = y0 - y1;
		double phi21 = x0 - x2;
		double phi22 = x1 - x0;

		
		// Assemblaggio matrice di stiffness
		Tmp(n*i0+i0) += (0.5)*invdetJ*((phi01*phi01) + (phi02*phi02));
		Tmp(n*i0+i1) += (0.5)*invdetJ*((phi11*phi01) + (phi12*phi02));
		Tmp(n*i0+i2) += (0.5)*invdetJ*((phi21*phi01) + (phi22*phi02));		
		Tmp(n*i1+i0) += (0.5)*invdetJ*((phi01*phi11) + (phi02*phi12));
		Tmp(n*i1+i1) += (0.5)*invdetJ*((phi11*phi11) + (phi12*phi12));		
		Tmp(n*i1+i2) += (0.5)*invdetJ*((phi21*phi11) + (phi22*phi12));	
		Tmp(n*i2+i0) += (0.5)*invdetJ*((phi01*phi21) + (phi02*phi22));	
		Tmp(n*i2+i1) += (0.5)*invdetJ*((phi11*phi21) + (phi12*phi22));	
		Tmp(n*i2+i2) += (0.5)*invdetJ*((phi21*phi21) + (phi22*phi22));
	}
	int j = 0;
	int k = 0;
	for(int i = 0; i < n*n; ++i) {
		if(Tmp(i) != 0.) {
			// Valore
			A.val(j) = Tmp(i);
			// Indice colonna
			A.col_ind(j) = i % n;
			// Indice riga
			if (i % n == 0) {
				A.row_ptr(k) = j;
				k++;
			}
			j++;
		}
	}
}
/*
 * Costruzione termine noto
 */
void Poisson::makeTermineNoto(Vector &F) {
	// Calcolo dimensione del vettore termine noto
	int n = 0;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	while (itV != cdt_.finite_vertices_end()) {
		if (!(itV->info().isDirichlet()))
			++n;
		++itV;
	}
	//Assegnazione del termine noto
	F.newsize(n);
	F = 1;
}


#endif // POISSON_HPP
