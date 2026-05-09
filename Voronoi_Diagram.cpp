#include "Voronoi_Diagram.h"
#if TEST==1
#include "doctest.h"
#endif
/**
   p is vector of generator points/ site of the voronoi diagram 
**/


CVoronoiDiagram::CVoronoiDiagram(const std::vector<Point> &p):
   mBbox(Point(),0),
   mGenPoints(p)
{
   int index = 0;
   std::vector<double> xVals;
	std::vector<double> yVals;
   for (const Point& i: mGenPoints) 
   {
      mEventQueue.insert(Event_s(i, index)); // construct the new event prio queue out of generator points vectors
      mFaces.push_back(Face_s(i)); // generate faces vectors out of generator points vectors
	  xVals.push_back(i.x);
	  yVals.push_back(i.y);
      index++;
   }
   // find max x, min x , max y and min y of all gen points
	double r = 150;  // distance from the boundary box to the actual minimal and maximal values
	double xMax = *std::max_element(xVals.begin(), xVals.end());
	double xMin = *std::min_element(xVals.begin(), xVals.end());
	double yMax = *std::max_element(yVals.begin(), yVals.end());
	double yMin = *std::min_element(yVals.begin(), yVals.end());
	// take the biggest distance as length of the boundary square 
	double distance = std::max(xMax - xMin, yMax - yMin) + 2*r;
	
	BoundaryBox_s bBox( Point(xMin - r, yMax + r ), distance );  // create boundary box
    mBbox = bBox;  // copy constructor
}

/*
 return array of generator points for drawing
*/
std::vector<Point> CVoronoiDiagram::getGenPoints()
{
	return  mGenPoints;
}

/*
 return array of edges
*/
std::vector<Half_edge_s> CVoronoiDiagram::getEdges()
{
	return mHalfEdges;
}

CVoronoiDiagram::BoundaryBox_s CVoronoiDiagram::getbBox()
{
    return mBbox;	 
}

/*
  returns highest prio event - largest  y
*/
Event_s CVoronoiDiagram::getTopPrioEvent()
{
   Event_s res = *mEventQueue.begin();
   mEventQueue.erase(mEventQueue.begin());
   return res;
}

void CVoronoiDiagram::insertNewEvent(Event_s event)
{
   mEventQueue.insert(event);
}

bool CVoronoiDiagram::isEventQueueEmpty()
{
   return mEventQueue.empty();
}


/**
 calculate areas of voronoi cells
 
 return: vector of areas, index is same as vector of faces.
 if area is infinite, -1 is returned.
 
       *  genPoint1
	   Origin---->
	   <----Dest
	   *  genPoint2
*/
std::vector<double> CVoronoiDiagram::calculateAreas()
{
	std::vector<double> resultVector;
   // iterate over faces vector
   int faceIndex = 0;
   double area = 0;
   int firstInd = -1;
   int currIndex = 0;
   for(auto face: mFaces)
   {
      // save first outer edge
	  area = 0;
      firstInd = face.outerComponentIndex;	  
      currIndex = firstInd;
	  do
	  {
	     // calcuate the base length of the triangle
         Point Origin = getEdge(currIndex).Origin;
		 Point Dest = getEdge(getTwinIndex(currIndex)).Origin;
		 double triangleBaseLen = vectorOps::euclidDist(Origin, Dest);
		 // calculate the height of the triangle
		 Point genPoint1 =  mGenPoints.at(getEdge(currIndex).pointIndex);
		 Point genPoint2 =  mGenPoints.at(getEdge(getTwinIndex(currIndex)).pointIndex);
		 double height =  vectorOps::euclidDist(genPoint1, genPoint2)/2;
		 area +=  0.5 * (height * triangleBaseLen);   // update area
		 // go to the next edge index
		 currIndex =  getEdge(currIndex).nextHalfEdgeIndex;
	  } while(currIndex != -1 || currIndex != firstInd); // as long  as there is a next edge
	  
	  if(currIndex == firstInd) // the area is finite
	  {
	     resultVector.push_back(area); 
	  } else { // edge list is not cyclical -> infinite area
	     resultVector.push_back(-1);
	  }
	  faceIndex++;
   }
   return resultVector;
}

/**
   extend the edges by th intersections with  the boundary box
*/
void CVoronoiDiagram::intersectEdgesWithBbox()
{
        int currEdgeIndex = 0; // needed to get twin edges etc.
	for(Half_edge_s& edge : mHalfEdges)
	{
	   // calculate intersections of voronoi line segments or rays with the boundary box
	   if(edge.nextHalfEdgeIndex != -1 && edge.prevHalfEdgeIndex != -1) // if it is a line segment, origin and destination are known
	   {
		  //Point orig = edge.Origin;
		  //Point dest =  mHalfEdges[edge.TwinEdge].Origin;
		  // do nothing, both ends of an edge are already set
	   } else {
		  Point genPoint1 = mGenPoints.at(edge.pointIndex);
		  int twinPointIndex = mHalfEdges.at(edge.TwinEdge).pointIndex;
		  Point genPoint2 =  mGenPoints.at(twinPointIndex);
		  Vector orthoVec(genPoint1.x - genPoint2.x, genPoint1.y - genPoint2.y);
		  Vector lineVec(-orthoVec.y, orthoVec.x);  // vector that is orthogonal to the vector that connects the generation points (genPoint2 to the left)
		  
	      if(edge.nextHalfEdgeIndex == -1 && edge.prevHalfEdgeIndex == -1)
	      {
			  Point origPoint( 0.5*(genPoint1.x + genPoint2.x), 0.5*(genPoint1.y + genPoint2.y) );  // this is the starting point for our vector
			  std::pair<Point, Point> points = mBbox.calcIntersectionsWithLine(origPoint, lineVec);  //  if inifinite ray in both directions, calculate two intersection points with the boundary box
			  Point orig = points.first;
		      Point dest =  points.second;
		      // set both points so we have a closed dcel and it can be drawed
		      edge.Origin = orig;
		      getEdge(getTwinIndex(currEdgeIndex)).Origin = dest;
	      } else {
			  // if inifinite ray in one direction, calculate one intersection  with boundary box
			  if(edge.nextHalfEdgeIndex != -1)
			  {
				  Point orig = mHalfEdges[edge.TwinEdge].Origin; // this is the starting point for our vector
				  std::pair<Point, Point> points = mBbox.calcIntersectionsWithLine(orig, lineVec);
				  Point dest = points.first; // take the destination point for negative param t (since the vector lineVec is in other direction)
				  //  set origin
				  edge.Origin = dest;
			  } else {
				  Point orig = edge.Origin; // this is the starting point for our vector
				  std::pair<Point, Point> points = mBbox.calcIntersectionsWithLine(orig, lineVec);
				  Point dest = points.second; // take the destination point for positive param t
				  // set destination
				  getEdge(getTwinIndex(currEdgeIndex)).Origin = dest;
			  }
		  }
	   }
	   currEdgeIndex++;
	}
	
}


/*
  create two new edges between twon points : push to edges vector and connect them with each other
  
  int edgePointindex, int twinEdgePointindex: indeces of faces/points in the face vector between which the edge pair is created
*/
void CVoronoiDiagram::createTwinEdges(int edgePointindex, int twinEdgePointindex)
{
   Half_edge_s newEdge(mFaces.at(edgePointindex), edgePointindex);
   Half_edge_s newEdgeTwin(mFaces.at(twinEdgePointindex), twinEdgePointindex);
   
   mHalfEdges.push_back(newEdge);
   mHalfEdges.push_back(newEdgeTwin);
   
   int newEdgeIndex =  mHalfEdges.size() - 2;
   int newEdgeTwinIndex =  mHalfEdges.size() - 1;
   
   mHalfEdges[newEdgeIndex].TwinEdge = newEdgeTwinIndex; // set twin edges 
   mHalfEdges[newEdgeTwinIndex].TwinEdge = newEdgeIndex;   
}

/*
  get edge at given index
*/
Half_edge_s& CVoronoiDiagram::getEdge(int edgeIndex)
{
   return mHalfEdges.at(edgeIndex);
}

/*
 get index of twin edge
*/
int CVoronoiDiagram::getTwinIndex(int edgeIndex)
{
   return mHalfEdges.at(edgeIndex).TwinEdge;
}


/*
  connect prevEdge <--> Edge
*/
 void CVoronoiDiagram::connectEdges(int prevEdge, int Edge)
 {
    getEdge(prevEdge).nextHalfEdgeIndex = Edge;
    getEdge(Edge).prevHalfEdgeIndex = prevEdge;
 }
 
 /**
 handle new vertex during circle event 
 Point vertex: new vertex
 int middleArcPointIndex: index of the site point of the deleted arc
 Arc_s* prevArc, Arc_s* succArc: predecessor and successor arcs of the deleted arc
 
 */
 void CVoronoiDiagram::handleNewVertex(Point vertex, int middleArcPointIndex, Arc_s* prevArc, Arc_s* succArc)
 {
   // 2. Add the center of the circle causing the event as a vertex record to the doubly-connected edge list D storing the Voronoi diagram under construc
   //-tion. Create two half-edge records corresponding to the new breakpoint of the beach line. Set the pointers between them appropriately. Attach the
   //three new records to the half-edge records that end at the vertex.
   
   // we assume here that arcs trace out edges that belong to their own sitePoints
   getEdge(prevArc->rightEdgeIndex).Origin = vertex; // set vertex point that connects the new edges: prevArc <--rightedge--- vertex ---leftedge--> succArc
   getEdge(getTwinIndex(succArc->leftEdgeIndex)).Origin = vertex;  // set the origin of the twin of the successor edge
   connectEdges( getTwinIndex(prevArc->rightEdgeIndex), getTwinIndex(succArc->leftEdgeIndex) ); // connect rightTwin to leftTwin, these are edges of the disappearing arcs face
   // new breakpoint of the beach line is between prev and next arc -> new edge
   createTwinEdges(prevArc->pointIndex,  succArc->pointIndex);  // Create two half-edge records corresponding to the new breakpoint of the beach line
   
   int newEdgeIndex =  mHalfEdges.size() - 2;  // prev Arc`s  face
   int newEdgeTwinIndex =  mHalfEdges.size() - 1;
   // connect newEdge and newEdgeTwin to the edges of the face of disappearing arc      
   
   connectEdges(newEdgeIndex, prevArc->rightEdgeIndex);
   connectEdges(succArc->leftEdgeIndex, newEdgeTwinIndex);
   
   //Update faces with outer component edges
   mFaces.at(prevArc->pointIndex).outerComponentIndex = newEdgeIndex;
   mFaces.at(succArc->pointIndex).outerComponentIndex = newEdgeTwinIndex;
   mFaces.at(middleArcPointIndex).outerComponentIndex = getTwinIndex(succArc->leftEdgeIndex);   // face of the site in the middle
   getEdge(newEdgeTwinIndex).Origin  = vertex;  // set edge connections in the dcel   : prevArc <--rightedge--- vertex ---leftedge--> succArc
                                                                                                  //              |   
																							     //             newEdge
   succArc->leftEdgeIndex = newEdgeIndex;  // save new currently traced out twin edges , the "old" edges should be accessible from their faces and the vector of edges
   prevArc->rightEdgeIndex = newEdgeIndex;
 }
 
 void CVoronoiDiagram::checkTriplesForConvergenceWithEvolution(Arc_s* leftMiddle, Arc_s* rightMiddle, double sweepLineY)
 {
   Point lowestCirclePointLeft; 
   Point lowestCirclePointRight; 
   Point intersectionPointLeft;
   Point intersectionPointRight;
   
   bool leftArcConvergence = mBeachBST.doCurrentBreakPointsConverge(leftMiddle->prev, leftMiddle, leftMiddle->next, lowestCirclePointLeft, intersectionPointLeft, sweepLineY );
   bool rightArcConvergence = mBeachBST.doCurrentBreakPointsConverge(rightMiddle->prev, rightMiddle, rightMiddle->next, lowestCirclePointRight, intersectionPointRight, sweepLineY );
   
   // for a circle event the event point that we store is the lowest point of the circle, with a pointer to the leaf in T that represents the arc that will disappear in the event.
   
   //add pointers between the node in T and the node in Q. Do the same for the
   //triple where the new arc is the right arc
   if(leftArcConvergence)
   {
      Event_s event(lowestCirclePointLeft, -1);  // circle event -> no such index exists in the vector of gen points
      event.isSiteEvent = false; // true->site event, false -> circle event
      event.circleEventArc = leftMiddle;
      event.intersectionPoint  = intersectionPointLeft;
      leftMiddle->isEventSet = true;
      // also set intersection point
      leftMiddle->circleEvent =  event;
      insertNewEvent(event);
   }
   
   if(rightArcConvergence)
   {
      Event_s event(lowestCirclePointRight, -1);
      event.isSiteEvent = false;
      event.circleEventArc = rightMiddle;
      event.intersectionPoint  = intersectionPointRight;
      rightMiddle->isEventSet = true;
      rightMiddle->circleEvent =  event;
      insertNewEvent(event);
   } 
	 
 }

void CVoronoiDiagram::HandleSiteEvent(Event_s event)
{
   // else look for arc above the new site point
   if( mBeachBST.isEmpty() )
   {
      // if BST is empty, insert new arc based on new site point into the binary search tree
      Arc_s* arc = new Arc_s(event.point, event.pointIndex);
      mBeachBST.insertArcAfter(arc, nullptr);
      return;
   }
   Arc_s* arc = new Arc_s(event.point, event.pointIndex);
   double sweepLineYpos = event.point.y;
   Arc_s* arcAbove = mBeachBST.findArcAboveInLinearTime(arc, sweepLineYpos);  //TODO special handling needed if no arc above exists!!

   if(arcAbove->isEventSet) // check pointer to a circle event from arcAbove to  the priority queue and delete it 
   {
      arcAbove->isEventSet = false;
      auto erase_cnt = mEventQueue.erase( arcAbove-> circleEvent ); // delete circle event, O(log n) in size n of event queue complexity
      // check erase count is 1
      assert(erase_cnt == 1);
   } 
   mBeachBST.insertArcAfter (arc, arcAbove);
   Arc_s* arcAboveCopy = new Arc_s(*arcAbove);  //               arcAbove--arc--arcAboveCopy
   // insert a copy  of the arc above as well 
   mBeachBST.insertArcAfter (arcAboveCopy, arc);  
   // Create new half-edge records in the Voronoi diagram structure for the edge separating V(pi ) and V(p j ), which will be traced out by the two new breakpoints
   createTwinEdges(event.pointIndex, arcAbove->pointIndex);
   int newEdgeIndex =  mHalfEdges.size() - 2;
   int newEdgeTwinIndex =  mHalfEdges.size() - 1;
   
   //Update faces with outer component edges
   mFaces.at(event.pointIndex).outerComponentIndex = newEdgeIndex;
   mFaces.at(arcAbove->pointIndex).outerComponentIndex = newEdgeTwinIndex;
   
   arcAboveCopy->leftEdgeIndex = newEdgeTwinIndex;  // every arc traces out the half edges tht belong to its face
   arcAbove->rightEdgeIndex = newEdgeTwinIndex; 
   
   arc->leftEdgeIndex = newEdgeIndex;  // new arc is in the middle and both breakpoints trace out the same edges, alway the half edge with new point = site point
   arc->rightEdgeIndex = newEdgeIndex;
  
   // Check the triple of consecutive arcs where the new arc for pi is the left arc to see if the breakpoints converge. If so, insert the circle event into Q and
   //checkTriplesForConvergence(arc->next, arc->prev);
   checkTriplesForConvergenceWithEvolution(arc->next, arc->prev, sweepLineYpos);
}
 
void CVoronoiDiagram::HandleCircleEvent(Event_s event)
{
   double sweepLineYpos = event.point.y;
   Arc_s*  arcToDisappear = event.circleEventArc;
   // before the arc is deleted,   Delete all circle events involving α from Q; these can be found using the pointers from the predecessor and the successor of γ in T
   Arc_s* prevArc = arcToDisappear->prev;
   Arc_s* succArc = arcToDisappear->next;
   if(prevArc->isEventSet)
   {
      prevArc->isEventSet = false;
      mEventQueue.erase( prevArc-> circleEvent );
   }
   if(succArc->isEventSet)
   {
      succArc->isEventSet = false;
      mEventQueue.erase( succArc-> circleEvent );  
   }
   
   int middleArcPointIndex = arcToDisappear->pointIndex; //TODO
   mBeachBST.deleteArc (arcToDisappear);
   delete arcToDisappear;
   //2. Add the center of the circle causing the event as a vertex record to the doubly-connected edge list D storing the Voronoi diagram under construc
    //-tion. Create two half-edge records corresponding to the new breakpoint of the beach line. Set the pointers between them appropriately. Attach the
    //three new records to the half-edge records that end at the vertex.
   Point vertex = event.intersectionPoint;
   handleNewVertex(vertex, middleArcPointIndex, prevArc, succArc);
   
   //3. Check the new triple of consecutive arcs that has the former left neighbor of α as its middle arc to see if the two breakpoints of the triple converge.
   //If so, insert the corresponding circle event into Q. and set pointers between the new circle event in Q and the corresponding leaf of T. Do the same for
   //the triple where the former right neighbor is the middle arc.
  
   //checkTriplesForConvergence(prevArc, succArc);
   
   checkTriplesForConvergenceWithEvolution(prevArc, succArc, sweepLineYpos);
   
}

/*Tests*/
#if TEST==1

TEST_CASE("CVoronoiDiagram_Constructor_Test") 
{
    std::vector<Point> p = { {0,1}, {2,3} };
    CVoronoiDiagram diagram(p);
    CHECK(diagram.isEventQueueEmpty() == false );
}

TEST_CASE("CVoronoiDiagram_PriorityQueueLargestYFirstCheck") 
{
   std::vector<Point> p = { {0,1}, {2,3}, {-3, -2} };
   CVoronoiDiagram diagram(p);
   CHECK(diagram.getTopPrioEvent().point.y == 3 );
}

TEST_CASE("CVoronoiDiagram_PriorityQueueDeleteValCheck") 
{
   std::vector<Point> p = { {0,1}, {2,3}, {-3, -2} };
   CVoronoiDiagramTestInterface diagram(p);
   Event_s event(Point(2,3), -1);
   CHECK(diagram.getTopPrioEvent().point.y == 3 );
   diagram.deleteEvent(event);
   CHECK(diagram.getTopPrioEvent().point.y == 1 );
}

TEST_CASE("CVoronoiDiagram_PriorityQueueElementRemovalTest") 
{
   std::vector<Point> p = { {0,1}, {2,3}, {-3, -2} };
   CVoronoiDiagram diagram(p);
   CHECK(3 ==  diagram.getEventQueueSize() );
   diagram.getTopPrioEvent();
   CHECK(2 ==  diagram.getEventQueueSize() );  
}

TEST_CASE("TwinEdgeCreationTest")
{
   std::vector<Point> p = { {0,1}, {2,3}, {-3, -2} };
   CVoronoiDiagramTestInterface diagram(p);
   diagram.createTwinEdgesInt(0, 1);  // create edge pair between  {0,1} and {2,3}
   int twinEdgeInd = diagram.getEdgeVectorSize() - 1 ;
   int EdgeInd = diagram.getEdgeVectorSize() - 2 ;
   
   Half_edge_s& edge = diagram.getEdgeInt(EdgeInd); 
   Half_edge_s& edgeTwin = diagram.getEdgeInt(twinEdgeInd); 
   CHECK(edgeTwin.TwinEdge == EdgeInd);
   CHECK(edge.TwinEdge == twinEdgeInd);
}

TEST_CASE("EdgeValueSettingWithgetEdgeTest")
{
   std::vector<Point> p = { {0,1}, {2,3}, {-3, -2} };
   CVoronoiDiagramTestInterface diagram(p);
   diagram.createTwinEdgesInt(0, 1);  // create edge pair between  {0,1} and {2,3}
   int twinEdgeInd = diagram.getEdgeVectorSize() - 1 ;
   //int EdgeInd = diagram.getEdgeVectorSize() - 2 ;
   
   //Half_edge_s& edge = diagram.getEdgeInt(EdgeInd); 
   Half_edge_s& edgeTwin = diagram.getEdgeInt(twinEdgeInd);
   edgeTwin.TwinEdge = -1;
   CHECK(diagram.getEdgeInt(twinEdgeInd).TwinEdge == -1);
}


// arc1 <-> deleted_arc <-> arc2

/*
                                    p_j         
			  Edge            Edge
                         ----->   vertex  ---->  arc2
	               arc1 <----          <-----
			            ^  |
			            |  |
			
			     p_i    |  |     p_k
				    |  |
				      v
*/

TEST_CASE("handleNewVertex_VertexTest")
{
   std::vector<Point> p = { {1,0}, {2,3}, {3, -1} };  // the arc of the middle point is deleted
   CVoronoiDiagramTestInterface diagram(p); // init faces and site events
   Vertex_s intersectionPoint(2,5); // some values, not necessary plausible

   diagram.createTwinEdgesInt(0, 1);
   diagram.createTwinEdgesInt(1, 2);
   
   Arc_s arc1(p[0], 0);
   Arc_s arc2(p[2], 2);
   
   arc1.rightEdgeIndex = 0;  // arcs and the edges being traced out should have identical site points
   arc2.leftEdgeIndex = 3;  // thes edge indeces will be updated after handleNewVertex
   
   int rightEdgeIndex = 0;
   int leftEdgeIndex = 3;
  
   diagram.handleNewVertexInt(intersectionPoint, 1, &arc1, &arc2);
      
   CHECK(diagram.getEdgeInt(rightEdgeIndex).Origin.x == intersectionPoint.x); // origin of right edge of arc1 is now the new vertex
   CHECK(diagram.getEdgeInt(rightEdgeIndex).Origin.y == intersectionPoint.y);
   
   int leftEdgeIndexTwin =  diagram.getTwinIndexInt(leftEdgeIndex );
   CHECK(diagram.getEdgeInt(leftEdgeIndexTwin).Origin.x == intersectionPoint.x);   // destination of the left edge of arc2 is the new vertex
   CHECK(diagram.getEdgeInt(leftEdgeIndexTwin).Origin.y == intersectionPoint.y);
   	
   // new edge of the last point p_k should have vertex as origin
   int pk_new_edge = 5;
   CHECK(diagram.getEdgeInt(pk_new_edge).Origin.x == intersectionPoint.x); 
   CHECK(diagram.getEdgeInt(pk_new_edge).Origin.y == intersectionPoint.y); 
}

TEST_CASE("handleNewVertex_EdgeConnectionTest")
{
   std::vector<Point> p = { {1,0}, {2,3}, {3, -1} };  // the arc of the middle point is deleted
   CVoronoiDiagramTestInterface diagram(p); // init faces and site events
   Vertex_s intersectionPoint(2,5); // some values, not necessary plausible

   diagram.createTwinEdgesInt(0, 1);
   diagram.createTwinEdgesInt(1, 2);
   
   Arc_s arc1(p[0], 0);
   Arc_s arc2(p[2], 2);
   
   arc1.rightEdgeIndex = 0;  // arcs and the edges being traced out should have identical site points
   arc2.leftEdgeIndex = 3;
   
   int rightEdgeIndex = 0;
   int leftEdgeIndex = 3;
   diagram.handleNewVertexInt(intersectionPoint, 1, &arc1, &arc2);
	
   // prev edge of right edge of arc1 is the new edge of p_i
   int p_i_EdgeIndex = 4;
    //  check edge connections
   CHECK(diagram.getEdgesInt().size() == 6); // 2 new edges should be created
   CHECK(diagram.getEdgeInt(rightEdgeIndex).prevHalfEdgeIndex == p_i_EdgeIndex);
   // nexte edge of the new edge of p_i is right edge of arc1
   CHECK(diagram.getEdgeInt(p_i_EdgeIndex).nextHalfEdgeIndex == rightEdgeIndex);
   // next edge of twin of right edge of arc1 is  twin of the left edge of arc2
   int rightEdgeIndexTwin =  diagram.getTwinIndexInt(rightEdgeIndex);
   int leftEdgeIndexTwin =  diagram.getTwinIndexInt(leftEdgeIndex );
   CHECK(diagram.getEdgeInt(rightEdgeIndexTwin).nextHalfEdgeIndex  == leftEdgeIndexTwin);
   CHECK(diagram.getEdgeInt(leftEdgeIndexTwin).prevHalfEdgeIndex == rightEdgeIndexTwin);
   
   // next edge of the left edge of arc2 is new edge of p_k
   int p_k_EdgeIndex = 5;
   CHECK(diagram.getEdgeInt(leftEdgeIndex).nextHalfEdgeIndex  == p_k_EdgeIndex);
   CHECK(diagram.getEdgeInt(p_k_EdgeIndex).prevHalfEdgeIndex  == leftEdgeIndex);
}

TEST_CASE("LineSegmentIntersectionTest")
{
   CVoronoiDiagram::Line_segment_s lineSeg(Point(0,0), Point(1,1));
   auto res = lineSeg.calcIntersectionsWithLineVector(Point(1,0), Vector(-1,1));
   CHECK(res.first.x == 0.5); 
   CHECK(res.first.y == 0.5); 
}

TEST_CASE("BoundaryBoxIntersectionTest")
{
   CVoronoiDiagram::BoundaryBox_s bbox(Point(0, 1), 1 );
   auto res = bbox.calcIntersectionsWithLine(Point(0.5, 0.5), Vector(1,0));
   CHECK(res.first.x == 1); 
   CHECK(res.first.y == 0.5); 
   
   CHECK(res.second.x == 0); 
   CHECK(res.second.y == 0.5); 
}


TEST_CASE("checkTriplesForConvergenceEvolutionCircleEventTest")
{
   std::vector<Point> p = { Point()};  // the arc of the middle point is deleted
   CVoronoiDiagramTestInterface diagram(p); // init faces and site events
   diagram.getTopPrioEvent();
   CHECK(diagram.isEventQueueEmpty());
   Arc_s leftMiddlePrev (Point(1, 0), 0);  // triple of converging arcs
   Arc_s leftMiddle (Point(2, 1), 1);
   Arc_s leftMiddleNext( Point(3, 0), 2);
   
   leftMiddlePrev.next =  &leftMiddle;  // connect arcs
   leftMiddle.prev =  &leftMiddlePrev;
   
   leftMiddle.next = &leftMiddleNext;
   leftMiddleNext.prev = &leftMiddle;
   
   Arc_s right(Point(), 0);
   double sweeplineY = 0;
   
   diagram.checkTriplesEvolutionForConvergenceInt(&leftMiddle, &right, sweeplineY);
   
   REQUIRE (!diagram.isEventQueueEmpty());
   Event_s event = diagram.getTopPrioEvent();
   REQUIRE (!event.isSiteEvent ); // check if circle event
   
   REQUIRE ( approxEquality_n::isEqualZero(event.point.x - 2 ) ); // check lowest circle point
   REQUIRE ( approxEquality_n::isEqualZero(event.point.y  + 1) );
   
   REQUIRE ( approxEquality_n::isEqualZero(event.intersectionPoint.x - 2) ); // check  intersection point
   REQUIRE ( approxEquality_n::isEqualZero(event.intersectionPoint.y) );
}



// Integration Tests

TEST_CASE("IntegrationTest2Points_2HalfEdges")
{
   std::vector<Point> vector = { Point(100,200), Point(300,300) };
   CVoronoiDiagramTestInterface vorDiagram(vector);
   
   while(!vorDiagram.isEventQueueEmpty())
   {
      Event_s event = vorDiagram.getTopPrioEvent();
      if(event.isSiteEvent)
      {
         vorDiagram.HandleSiteEvent(event);
      } else {
         vorDiagram.HandleCircleEvent(event);
     }
   }
   std::vector<Half_edge_s> edges = vorDiagram.getEdgesInt();
   REQUIRE(edges.size() == 2); 
}


TEST_CASE("IntegrationTest3Points_CircleEventCorrect")
{
   std::vector<Point> vector = { Point(100,200), Point(200,300) , Point(300,205) };
   CVoronoiDiagramTestInterface vorDiagram(vector);
   while(!vorDiagram.isEventQueueEmpty())
   {
      Event_s event = vorDiagram.getTopPrioEvent();
      if(event.isSiteEvent)
      {
         vorDiagram.HandleSiteEvent(event);
      } else {
             REQUIRE(event.isSiteEvent == false); 
	     REQUIRE( vorDiagram.isEventQueueEmpty() ); 
	  }
   }
   
}

TEST_CASE("IntegrationTest3Points_6HalfEdges_1")
{
   std::vector<Point> vector = { Point(100,200), Point(200,300) , Point(300,205) };
   CVoronoiDiagramTestInterface vorDiagram(vector);
   while(!vorDiagram.isEventQueueEmpty())
   {
      Event_s event = vorDiagram.getTopPrioEvent();
      if(event.isSiteEvent)
      {
         vorDiagram.HandleSiteEvent(event);
      } else {
         vorDiagram.HandleCircleEvent(event);
     }
   }
   std::vector<Half_edge_s> edges = vorDiagram.getEdgesInt();
   CHECK(edges.size() == 6); 
}

TEST_CASE("IntegrationTest3Points_6HalfEdges_2")
{
   std::vector<Point> vector = { Point(160,220) , Point(220,205), Point(230,225) };
   CVoronoiDiagramTestInterface vorDiagram(vector);
   while(!vorDiagram.isEventQueueEmpty())
   {
      Event_s event = vorDiagram.getTopPrioEvent();
      if(event.isSiteEvent)
      {
         vorDiagram.HandleSiteEvent(event);
      } else {
         vorDiagram.HandleCircleEvent(event);
     }
   }
   std::vector<Half_edge_s> edges = vorDiagram.getEdgesInt();
   CHECK(edges.size() == 6); 
}

TEST_CASE("IntegrationTest4Points_1") // 2 verteces are  connected by an edge
{
   std::vector<Point> vector = { Point(160,220), Point(200,300) , Point(220,205), Point(230,225) };
   CVoronoiDiagramTestInterface vorDiagram(vector);
   while(!vorDiagram.isEventQueueEmpty())
   {
      Event_s event = vorDiagram.getTopPrioEvent();
      if(event.isSiteEvent)
      {
         vorDiagram.HandleSiteEvent(event);
      } else {
         vorDiagram.HandleCircleEvent(event);
     }
   }
   std::vector<Half_edge_s> edges = vorDiagram.getEdgesInt();
   CHECK(edges.size() == 10); 
}

TEST_CASE("IntegrationTest4Points_2")  // 2 verteces are not connected by an edge
{
   std::vector<Point> vector = { Point(160,220), Point(200,300) , Point(160,420), Point(200,500) };
   CVoronoiDiagramTestInterface vorDiagram(vector);
   while(!vorDiagram.isEventQueueEmpty())
   {
      Event_s event = vorDiagram.getTopPrioEvent();
      if(event.isSiteEvent)
      {
         vorDiagram.HandleSiteEvent(event);
      } else {
         vorDiagram.HandleCircleEvent(event);
     }
   }
   std::vector<Half_edge_s> edges = vorDiagram.getEdgesInt();
   CHECK(edges.size() == 10); 
}

#endif
