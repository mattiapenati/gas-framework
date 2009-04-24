/* BibTeX
@article{0501496v2,
	eprint = {math/0501496v2},
	author = {Mark A. Taylor, Beth A. Wingate, Len P. Bos},
	title = {Several new quadrature formulas for polynomial integration in the triangle},
	year = {2007}
}
*/

template<typename Geometry, unsigned int Degree>
struct NewtonCotes_2 : public Method_2 < NewtonCotes_2 < Geometry , Degree > , Geometry > {
	/* nodi e pesi */
	static const double x[];
	static const double y[];
	static const double w[];
	/* numbero di punti */
	static const unsigned int nPoints;
};

/* Ordine 1 su triangolo : grado di esattezza 2 */
template < >
const unsigned int NewtonCotes_2 < Geometry::Triangle , 1 >::nPoints = 3;

template < >
const double NewtonCotes_2 < Geometry::Triangle , 1 >:: x[3] = {
	0.1666666666667,
	0.6666666666667,
	0.1666666666667
};

template < >
const double NewtonCotes_2 < Geometry::Triangle , 1 >:: y[3] = {
	0.6666666666667,
	0.1666666666667,
	0.1666666666667
};

template < >
const double NewtonCotes_2 < Geometry::Triangle , 1 >:: w[3] = {
	0.6666666666667,
	0.6666666666667,
	0.6666666666667
};

/* Ordine 2 su triangolo : grado di esattezza 4 */
template < >
const unsigned int NewtonCotes_2 < Geometry::Triangle , 2 >::nPoints = 6;

template < >
const double NewtonCotes_2 < Geometry::Triangle , 2 >:: x[6] = {
	0.0915762135098,
	0.8168475729805,
	0.0915762135098,
	0.1081030181681,
	0.4459484909160,
	0.4459484909160
};

template < >
const double NewtonCotes_2 < Geometry::Triangle , 2 >:: y[6] = {
	0.0915762135098,
	0.0915762135098,
	0.8168475729805,
	0.4459484909160,
	0.1081030181681,
	0.4459484909160
};

template < >
const double NewtonCotes_2 < Geometry::Triangle , 2 >:: w[6] = {
	0.2199034873106,
	0.2199034873106,
	0.2199034873106,
	0.4467631793560,
	0.4467631793560,
	0.4467631793560
};
