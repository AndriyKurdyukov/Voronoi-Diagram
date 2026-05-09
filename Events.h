 #ifndef EVENTS_H
#define EVENTS_H

#include <set>
#include "dcel.h"

/**
events and arcs declarations
*/ 

struct Arc_s;  // forward declaration

  struct Event_s  // two event types
  {
    bool isSiteEvent = true;  // true if site event, false if circle event
    Point point; // maybe an index to vector of faces?
    Point intersectionPoint; //  circle intersection point during circle event
    int pointIndex = -1; // generator point index in the face or in the gen point array,    -1 ->not set-> circle event
    Event_s(Point  newpoint, int index):
      point(newpoint),
      intersectionPoint(),
      pointIndex(index)
    {}
	        
    bool operator==(const Event_s& rhs) const
    {
       return (point.y == rhs.point.y && point.x == rhs.point.x); 
    } 
    Arc_s* circleEventArc = nullptr; 
  };
  
  struct Cmp // prio compare for the event queue
  {
     bool operator()(const Event_s& rhs, const Event_s& lhs) const
     {
        return lhs.point.y < rhs.point.y;
     } 
  };
  
  enum class Color_e
   {
      BLACK_e = 1, 
      RED_e
   };
   struct Arc_s
   {
      // red-black tree relevant
      Arc_s* parent=nullptr;
      Arc_s* leftChild=nullptr;
      Arc_s* rigthChild=nullptr;
      Color_e color;
      //double linked list relevant
      Arc_s* next=nullptr;
      Arc_s* prev=nullptr;
			  
      //Half_edge_s* leftEdge = nullptr;  // edge being traced out by the left breakpoint of the arc
      //Half_edge_s* rightEdge = nullptr;	  
      
      Arc_s(const Point point, int index): 
	     sitePoint(point), 
         pointIndex(index),
         circleEvent(point, -1)  // init to any value
      {}
      Arc_s(const Arc_s& arc):   // copy constructor
        sitePoint(arc.sitePoint),
	    pointIndex(arc.pointIndex),
	    circleEvent(arc.sitePoint, -1), // init to any value
	    leftEdgeIndex (arc.leftEdgeIndex),
	    rightEdgeIndex (arc.rightEdgeIndex)
      {}
      
      // doubly connected edge list relevant
      const Point sitePoint;  // TODO   should not be a reference, point refernce in the constructor  has a limited lifetime
      const int pointIndex = -1; // generator point index in the face or in the gen point array,   -1 ->not set
      Event_s circleEvent; 
      //circle event  Event_s where this arc disappears
      bool isEventSet = false; // is the circle event set or not
      int leftEdgeIndex = -1;
      int rightEdgeIndex = -1;
      
   };
  
#endif
