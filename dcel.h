#ifndef dcel
#define dcel
#include <cfloat>
#include <cmath>

//necessary components to implement double connected edge list as decribed in Computational Geometry by Mark de Berg · Otfried Cheong Marc van //Kreveld · Mark Overmars

struct Point // generator points , given by the coding challenge
{
  double x,y;
  Point() : x(0.0), y(0.0) {}
  Point(double x, double y) : x(x), y(y) {}
};

typedef Point Vertex_s;

typedef Point Vector;
struct Matrix_2x2 
{
   Vector colVector[2];
   // copy construtor needed?
};

namespace vectorOps
{
inline Point addVector (Point a, Vector b)
{
   return Point(a.x + b.x, a.y + b.y);   
}
   
inline Vector mult(Vector a, double c)
{
   return Vector(a.x * c, a.y * c);   
}
   
inline double euclidDist(Point a, Point b)
{
   return  std::sqrt( (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) ) ;
}
   
inline double euclidDistNorm(Vector a)
{
   return  std::sqrt( a.x*a.x + a.y*a.y ) ;
}

/*
multiply 2x2 matrix by 2x0 vector
return: 2x0 vector
*/
inline Vector multMatrixVector(Matrix_2x2 A, Vector vector)
{
   double l,t;
   l = A.colVector[0].x * vector.x + A.colVector[1].x * vector.y;
   t = A.colVector[0].y * vector.x + A.colVector[1].y * vector.y;
   return Vector (l,t);
}

inline Matrix_2x2 multMatrixScalar(Matrix_2x2 A, double scalar)
{
	 Matrix_2x2 res {{ mult(A.colVector[0], scalar), mult(A.colVector[1], scalar)   }};
	 return res;
}

/*

calculate adjugate of a 2x2 matrix
a b          d   -b
c d    ->   -c    a
*/

inline Matrix_2x2  calcAdjMatrix(Matrix_2x2 mat)
{
	Matrix_2x2 res{ {  { mat.colVector[1].y,   -mat.colVector[0].y  }, { -mat.colVector[1].x,  mat.colVector[0].x  }  }   };
	return res;
}

/*
  calculate matrix determinant
*/
inline double calcMatrixDet(Matrix_2x2 mat)
{
	return mat.colVector[0].x * mat.colVector[1].y - mat.colVector[0].y * mat.colVector[1].x;
}

/*
  only valid if inverse exists -> det(mat) != 0
*/
inline Matrix_2x2 calcInverseMatrix(Matrix_2x2 mat)
{
	double det = calcMatrixDet(mat);
	Matrix_2x2 adjMat = calcAdjMatrix(mat);
	return multMatrixScalar(adjMat, 1/det);
}

};

namespace approxEquality_n
{
   double const M_Epsilon = 0.00001; 
   inline bool isEqualZero(double val)  // check if value is close enough to zero
   {
      return (std::abs(val) < M_Epsilon);
   }
};


struct Face_s;

struct Half_edge_s
{
  Vertex_s Origin; 
  int TwinEdge = -1;
  const Face_s* IncidentFace; // face to the LEFT of the half-edge, it never changes during the algorithm.   it is ok to use ptr because vector of faces does not change ToDO use index to the vector of faces instead
  const int pointIndex = -1;
  
  //Half_edge_s* NextEdge=nullptr;
  //Half_edge_s* PrevEdge=nullptr;
  int nextHalfEdgeIndex = -1;  // index in the vector of edges
  int prevHalfEdgeIndex = -1;
  
  Half_edge_s(const Face_s& face, const int index): 
	 Origin(DBL_MAX, DBL_MAX), // set origin to an unplausible value
	 IncidentFace(&face), 
	 pointIndex(index)
  {
  }
};


struct Face_s // we will have as many faces as generator points
{
  //Half_edge_s* OuterComponent=nullptr;
  //Half_edge_s* InnerComponent=nullptr;
  
  int outerComponentIndex = -1;  // index of half edge in the vector of half edges
  int innerComponentIndex = -1;
  
  const Point& GeneratorPoint;  // we want to know coordinates of the relevant voronoi site generator point for this face
  
  Face_s(const Point& genPoint): GeneratorPoint(genPoint)
  {}
  Face_s(int outerComponent, int innerComponent,const Point& generatorPoint): 
    outerComponentIndex(outerComponent), 
    innerComponentIndex(innerComponent), 
    GeneratorPoint(generatorPoint)
  {}
};




#endif
