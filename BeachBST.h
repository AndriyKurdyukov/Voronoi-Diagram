#ifndef BEACHBST
#define BEACHBST

#include "dcel.h"
#include "set"
#include <cmath>
#include "Events.h"
#include <cfloat>
#include <utility>
/*
this class implements a necessary binary seach tree data structure for the fortune`s algorithm beachline 
instead of breakpoint and leaves just arcs are used
The beachline is represented as a red-black tree of arcs.
Arcs are also connected as a doubly linked list to make traversing along the beachline possible
*/

class CBeachBST
{
public:
   CBeachBST();
   
   //TODO ~CBeachBST();   destructor should delete all arcs starting with the leftmost arc
   
	// insert arc into the tree and double linked list and call fixupInsert afterwards
   void insertArcAfter (Arc_s* arcToInsert, Arc_s* targetArc);
   void insertArcBefore (Arc_s* arcToInsert, Arc_s* targetArc); // analogous to insertArcAfter
   void deleteArc (Arc_s* arc);
   
   struct Breakpoint_s  // breakpoints of an arc
   {
      double leftBreakpoint;
      double rightBreakpoint;
   };
   
   std::pair<Vector, Vector> calcBrPntEvolutionDirection(Arc_s* arc, double startSweepLineYpos);
   
   double calculateParabolaY(Point p, double x, double sweepLine);
   Point calculateParabolaIntersection (Point leftSite, Point rightSite, double sweepLineYpos );
   Breakpoint_s calculateArcBreakpoints(Arc_s* arc, double sweepLineYpos);
   
   Arc_s* findArcAbove(Arc_s* newSite, double sweepLineYpos); // find arc that lies directly above the new site, the arc intersections depend on the sweep line y position
   Arc_s* findArcAboveInLinearTime(Arc_s* newSite, double sweepLineYpos); // for testing only
	
   bool isEmpty()  // is BST empty?
   {
      return (mRoot == nullptr);
   }
    
   bool doCurrentBreakPointsConverge(Arc_s* left, Arc_s* middle,Arc_s* right, Point& circleEventPoint , Point& intersectionPoint, double sweelplineY);
   
private:
  
   void fixupInsert(Arc_s* arc); // fixup after inserting a new arc so that all red-black properties of the tree are ensured
   void fixupDelete(Arc_s* arc);
   Arc_s* mRoot = nullptr;  // root of the red black tree	
   
   Arc_s* mLeftMostArc = nullptr;  // TODO only use for testing , unnecessary when red black tree is used later
};

#endif
