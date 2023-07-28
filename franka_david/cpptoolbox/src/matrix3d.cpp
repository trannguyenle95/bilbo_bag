#include "matrix.hpp"


namespace math{

	Matrix3D::Matrix3D() : _nArr(1)
	{
		_m3d = new Matrix[_nArr];
		_m3d = _m3d - 1;
		_m3d[1] = Matrix();
	}
	Matrix3D::Matrix3D(int nArr) : _nArr(nArr)
	{
		if (nArr < 1)
			throw Exception("Matrix3D: nArr must be greater than 0");
		_m3d = new Matrix[nArr];
		_m3d = _m3d - 1;
		for (int i = 1; i <= nArr; i++)
			_m3d[i] = Matrix();
	}
	Matrix3D::Matrix3D(int nArr, int nRow, Real val) : _nArr(nArr)
	{
		if (nArr < 1 || nRow < 1)
			throw Exception("Matrix3D: nArr, nRow, nCol must be greater than 0");
		_m3d = new Matrix[nArr];
		_m3d = _m3d - 1;
		for (int i = 1; i <= nArr; i++)
			_m3d[i] = Matrix(nRow, nRow, val);
	}
	Matrix3D::Matrix3D(int nArr, int nRow, int nCol, Real val) : _nArr(nArr)
	{
		if (nArr < 1 || nRow < 1 || nCol < 1)
			throw Exception("Matrix3D: nArr, nRow, nCol must be greater than 0");
		_m3d = new Matrix[nArr];
		_m3d = _m3d - 1;
		for ( int i = 1; i <= nArr; i++)
			_m3d[i] = Matrix(nRow, nCol, val);
	}
	Matrix3D::Matrix3D(int nArr, int nRow, int nCol, const string &str) : _nArr(nArr)
	{
		if (nArr < 1 || nRow < 1 || nCol < 1)
			throw Exception("Matrix3D: nArr, nRow, nCol must be greater than 0");
		_m3d = new Matrix[nArr];
		_m3d = _m3d - 1;
		for ( int i = 1; i <= nArr; i++)
			_m3d[i] = Matrix(nRow, nCol, str);
	}
	Matrix3D::Matrix3D(const Matrix3D & orig)
	{
		copy(orig);
	}
	Matrix3D::~Matrix3D()
	{
		if (_m3d){
			_m3d = _m3d + 1;
			delete[] _m3d;
		}
	}
	void Matrix3D::copy(const Matrix3D & orig)
	{
		_nArr = orig.getnArr();
		_m3d = new Matrix[_nArr];
		_m3d = _m3d - 1;
		for (int i = 1; i <= _nArr; i++)
			_m3d[i] = orig(i);
	}
	Matrix3D & Matrix3D::operator=(const Matrix3D & orig)
	{
		if (_m3d){
			_m3d = _m3d + 1;
			delete[] _m3d;
		}
		copy(orig);
		return *this;
	}
	Matrix3D & Matrix3D::operator=(Real val)
	{
		if (!_m3d)
			throw Exception("Matrix3D::operator=(Real val) : Invalid Matrix");
		for (int i = 1; i <= _nArr; i++)
			_m3d[i] = val;
		return *this;
	}
	const Matrix & Matrix3D::operator()(int i)const
	{
		if (i < 1 || i > _nArr)
			throw Exception("Matrix3D: nArr must be greater than 0 and less than nArr");
		return _m3d[i];
	}
	Matrix & Matrix3D::operator()(int i)
	{
		if (i < 1 || i > _nArr)
			throw Exception("Matrix3D: nArr must be greater than 0 and less than nArr");
		return _m3d[i];
	}
	const Real & Matrix3D::operator()(int i, int j, int k)const
	{
		if (i < 1 || i > _nArr)
			throw Exception("Matrix3D: nArr must be greater than 0 and less than nArr");
		return _m3d[i](j,k);
	}
	Real & Matrix3D::operator()(int i, int j, int k)
	{
		if (i < 1 || i > _nArr)
			throw Exception("Matrix3D: nArr must be greater than 0 and less than nArr");
		return _m3d[i](j,k);
	}

}
