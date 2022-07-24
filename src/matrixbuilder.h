#ifndef MATRIX_BUILDER_H
#define MATRIX_BUILDER_H
#include "gauss.h"
#include "includes.h"
using namespace Wedge;

struct ConstantMatrixBuilder {
	int r,c;
	ex value;
	int rows() const {return r;}
	int cols() const {return c;}
	ex at(int i,int j) const {return value;}
};

template<typename Builder>
struct TransposeMatrixBuilder {
	const Builder& m;
	int rows() const {return m.cols();}
	int cols() const {return m.rows();}
	ex at(int i,int j) const {return m.at(j,i);}
};

template<typename Builder>
auto transpose(const Builder& m) {
	return TransposeMatrixBuilder<Builder>{m};
}

template<typename BuilderA,typename BuilderB>
struct ProductBuilder {
	const BuilderA& A;
	const BuilderB& B;
	int rows() const {return A.rows();}
	int cols() const {return B.cols();}
	ex at(int i,int j) const {
		ex sum;
		for (int h=0;h<A.cols();++h) sum+=A.at(i,h)*B.at(h,j);
		return sum;
	}
};

template<typename BuilderA,typename BuilderB> 
auto matrix_product(const BuilderA& A,	const BuilderB& B) {
	return ProductBuilder<BuilderA,BuilderB>{A,B};
}

struct ColumnVectorMatrixBuilder {
	const exvector& x;
	int rows() const {return x.size();}
	int cols() const {return 1;}
	ex at(int i,int j) const {return x[i];}
};

struct MatrixBuilder {
	const matrix& m;
	int rows() const {return m.rows();}
	int cols() const {return m.cols();}
	ex at(int i,int j) const {return m(i,j);}
};

template<typename Unknown>
struct UnknownsMatrixBuilder{
	exvector unknowns;
public:
	UnknownsMatrixBuilder(Name& name, int n) : unknowns{generate_variables<Unknown>(name,n)} {}
	int rows() const {return unknowns.size();}
	int cols() const {return 1;}
	ex at(int i,int j) const {
		assert(j==0);
		return unknowns[i];
	}
};

template<typename MatrixBuilder>
pair<int,int> adjoined_matrix_size(MatrixBuilder builder) {
	return make_pair(builder.rows(),builder.cols());
}

template<typename MatrixBuilder, typename... MatrixBuilders>
pair<int,int> adjoined_matrix_size(MatrixBuilder builder, MatrixBuilders... builders) {
	auto size=adjoined_matrix_size(builders...);
	if (builder.rows()!=size.first) throw std::invalid_argument("cannot adjoin matrices with different number of rows");
	size.second+= builder.cols();
	return size;
}

template<typename Matrix, typename MatrixBuilder>
Matrix fill(Matrix&& matrix, int from_column, MatrixBuilder builder) {
	for (int i=0;i<builder.rows();++i)
	for (int j=0;j<builder.cols();++j)
		matrix(i,j+from_column)=builder.at(i,j);
	return move(matrix);
}

template<typename Matrix, typename MatrixBuilder, typename... MatrixBuilders>
Matrix fill(Matrix&& matrix, int from_column, MatrixBuilder builder, MatrixBuilders... builders) {
	return fill(
		fill(move(matrix),from_column+builder.cols(),builders...), from_column,builder);
}

template<typename... MatrixBuilders>
Matrix adjoin(MatrixBuilders... builders) {
	auto size=adjoined_matrix_size(builders...);
	return fill(Matrix{size.first,size.second},0,builders...);
}

template<typename MatrixBuilder>
Matrix to_matrix(MatrixBuilder matrix_builder) {
	return fill(Matrix{matrix_builder.rows(),matrix_builder.cols()},0,matrix_builder);
}
#endif
