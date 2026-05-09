#ifndef VORONOI_DIAGRAM
#define VORONOI_DIAGRAM

#include <utility>
#include "dcel.h"
#include "vector"
#include "set"
#include "BeachBST.h"
#include "Events.h"
#include <algorithm>
#include <array>
#include <iostream>
#include <SFML/Graphics.hpp>
#include <cassert>

//necessary components to implement voronoi diagram realisation as decribed in Computational Geometry by Mark de Berg · Otfried Cheong Marc van Kreveld · Mark Overmars


class CVoronoiDiagram
{
  public:
  CVoronoiDiagram(const std::vector<Point> &p);
    
  Event_s getTopPrioEvent(); // todo get highest priority event from Q
  /**
    functions as described in computational geometry book
  */
  void HandleSiteEvent(Event_s event);
  void HandleCircleEvent(Event_s event);
  bool isEventQueueEmpty();
  int getEventQueueSize()
  {
     return mEventQueue.size();
  }
  
  std::vector<double> calculateAreas();
  
  void intersectEdgesWithBbox();
  
  struct Line_segment_s
  {
    Point mBegin;
    Point mEnd;
    Line_segment_s(Point begin, Point end):
       mBegin(begin),
       mEnd(end)
    {
    }
    
    /*
      return: Intersection point with the line segment, vector parameter with which the linevec is multiplied to produce an inetrsection
      if no point exists ( [DBL_MAX, DBL_MAX], DBL_MAX ) is returned
    */
    std::pair<Point, double> calcIntersectionsWithLineVector(Point startingPoint, Vector lineVec)
    {
      Vector lineSegmentVec(mEnd.x - mBegin.x, mEnd.y - mBegin.y );
      Matrix_2x2 A{{  {lineVec.x, lineVec.y },  {lineSegmentVec.x, lineSegmentVec.y }   }};
      Vector OrigPointDiff(mBegin.x - startingPoint.x, mBegin.y - startingPoint.y);
       // calculate inverse A times vector OrigPointDiff -> equals (l , -t)
	  Matrix_2x2 InvA;
	  //std::cout <<"Matrix det: "<<vectorOps::calcMatrixDet(A) <<"\n";
	  if(vectorOps::calcMatrixDet(A) != 0)
	  {
	     InvA = vectorOps::calcInverseMatrix(A);
	  } else {
		 return std::make_pair( Point(DBL_MAX,DBL_MAX), DBL_MAX ); ; //  no intersection since inverse does not exist
	  }
      Vector l_t = vectorOps::multMatrixVector(InvA, OrigPointDiff);  // x is parameterer of linevec l, y is param of line segment vector t
	  l_t.y = -l_t.y;
	  Point intersectionPoint;
	  if( l_t.y > 1 || l_t.y < 0 )  // is the intersection outside of the line segment 
	  {
	          //std::cout <<"l: "<<l_t.y <<"t: "<< l_t.x <<"\n";
		  // no intersection with the line segment exists
		  return std::make_pair( Point(DBL_MAX,DBL_MAX), DBL_MAX );
	  }
	  intersectionPoint = vectorOps::addVector(  startingPoint , vectorOps::mult( lineVec, l_t.x )) ;
	  return std::make_pair( intersectionPoint, l_t.x);
    }
  
  };
  
    /*
       mA - mB
	   |     |
	   mD - mC
  */
  struct BoundaryBox_s
  {
     double mLen = 0;	 
	 Point mA;
	 Point mB;
	 Point mC;
	 Point mD;
	 std::array<Line_segment_s, 4> mLineSegments; // line segments of the boundary box
	 
     BoundaryBox_s(Point leftUpperVertex, double len):
            mLen(len),
	    mA(leftUpperVertex),
		mB( mA.x + len , mA.y),
		mC( mA.x + len , mA.y - len),
	    mD( mA.x , mA.y - len),
		mLineSegments{ Line_segment_s(mA, mB), Line_segment_s(mB, mC), Line_segment_s(mC, mD), Line_segment_s(mD, mA)  }
	 {}
	 
	 /**
	   calculate intersection points with the boundary box of line-vector for  parameter t
	   return: {intersection point for t>= 0, intersection point for t < 0 }
	   [DBL_MAX, DBL_MAX] if no intersection exists
	 */
	 std::pair<Point, Point> calcIntersectionsWithLine(Point startingPoint, Vector lineVec)
	 {
		 auto result = std::make_pair(Point(DBL_MAX, DBL_MAX), Point(DBL_MAX, DBL_MAX)); 
		 for(Line_segment_s &lineSegment: mLineSegments)
		 {
			 auto res = lineSegment.calcIntersectionsWithLineVector(startingPoint, lineVec);
			 if(res.second >= 0 && res.second!= DBL_MAX)
			 {
				 result.first = res.first;  // intersection point for positive t
			 } else {
				 if(res.second < 0)
				 {
					 result.second = res.first; // intersection point for negative t
				 }
			 }
		 }
		 return result;
	 }
  };

  std::vector<Point> getGenPoints();
  std::vector<Half_edge_s> getEdges();
  BoundaryBox_s getbBox();
  
  protected:
  
  BoundaryBox_s mBbox;
  
  Half_edge_s& getEdge(int edgeIndex);
  int getTwinIndex(int edgeIndex);
  void createTwinEdges(int edgePointindex, int twinEdgePointindex);
  void connectEdges(int prevEdge, int Edge);
  void handleNewVertex(Point vertex, int middleArcPointIndex, Arc_s* prevArc, Arc_s* succArc);
  void checkTriplesForConvergenceWithEvolution(Arc_s* leftMiddle, Arc_s* rightMiddle, double sweepLineY);
  void insertNewEvent(Event_s event);
  std::set<Event_s, Cmp> mEventQueue;  // event queue  Q wth site and circle events, std::set to make deletion of obejcts possible
  std::vector<Face_s> mFaces;   // no new faces will be added after program start, but faces themselves will be modified
  const std::vector<Point>& mGenPoints;  // generator points
  std::vector<Half_edge_s> mHalfEdges; // doubly-connected edge list D
  CBeachBST mBeachBST;   // beach line represented by a balanced binary search tree T
};

#if TEST==1

class CVoronoiDiagramTestInterface: public CVoronoiDiagram   // interface class to access protected functions, needed for testing
{
   public:
   CVoronoiDiagramTestInterface(const std::vector<Point> &p):
	   CVoronoiDiagram(p)
	{}
	
	Half_edge_s& getEdgeInt(int edgeIndex)
	{
	   return getEdge(edgeIndex);
	}
	
	int getEdgeVectorSize() // get current number of elements i  the edge vector
	{
	   return mHalfEdges.size();
	}
	
	void createTwinEdgesInt(int edgePointindex, int twinEdgePointindex)
	{
	   createTwinEdges(edgePointindex, twinEdgePointindex);
	}
	
	void handleNewVertexInt(Point vertex, int middleArcPointIndex, Arc_s* prevArc, Arc_s* succArc)
	{
	   handleNewVertex(vertex, middleArcPointIndex, prevArc, succArc);
	}
	
	std::vector<Half_edge_s>  getEdgesInt()
	{
	   return  mHalfEdges;
	}
	
	
	int getTwinIndexInt(int edgeIndex)
	{
	   return getTwinIndex(edgeIndex);
	}
		
	void checkTriplesEvolutionForConvergenceInt(Arc_s* leftMiddle, Arc_s* rightMiddle, double sweeplineY)
	{
		checkTriplesForConvergenceWithEvolution(leftMiddle, rightMiddle,sweeplineY);
	}
	
	void deleteEvent(Event_s event)
	{
	   mEventQueue.erase( event );  
	}
	

		
	
};
#endif

#endif
