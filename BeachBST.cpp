#include "BeachBST.h"
#if TEST==1
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#endif

using namespace vectorOps;

CBeachBST::CBeachBST()
{}

/*
insert arc arcToInsert after   targetArc
*/
void CBeachBST::insertArcAfter (Arc_s* arcToInsert, Arc_s* targetArc)
{
   if(mRoot == nullptr)  // if tree is empty, create new root
   {
      mRoot = arcToInsert;
      mLeftMostArc = arcToInsert;
      arcToInsert->color = Color_e::BLACK_e; // root should always be black
      return;
   }
   // inserting into linked list
   if(targetArc->next != nullptr)
   {
      arcToInsert->next = targetArc->next;
      targetArc->next->prev = arcToInsert;
   }
   targetArc->next = arcToInsert;
   arcToInsert->prev = targetArc;
   // inserting into tree, not yet balanced
   arcToInsert->rigthChild = targetArc->rigthChild;
   targetArc->rigthChild = arcToInsert;
   arcToInsert->parent = targetArc;
   arcToInsert->color = Color_e::RED_e; // nodes to insert should be red
   fixupInsert(arcToInsert); // fix possible red-black property violations
}

void CBeachBST::insertArcBefore (Arc_s* arcToInsert, Arc_s* targetArc)
{
   if(mRoot == nullptr)  // if tree is empty, create new root
   {
      mRoot = arcToInsert;
	  mLeftMostArc = arcToInsert;
      arcToInsert->color = Color_e::BLACK_e;
      return;
   }
   
   // inserting into linked list
   arcToInsert->next = targetArc;
   targetArc->prev->next = arcToInsert;
   
   arcToInsert->prev = targetArc->prev;
   targetArc->prev = arcToInsert;
   
    // inserting into tree, not yet balanced
   arcToInsert->leftChild = targetArc->leftChild;
   targetArc->leftChild = arcToInsert;
   arcToInsert->parent = targetArc;
   arcToInsert->color = Color_e::RED_e; // nodes to insert should be red
   fixupInsert(arcToInsert); // fix possible red-black property violations
}

void CBeachBST::deleteArc (Arc_s* arc)
{
   //  delete node
   // unlink from the linked list
   if(arc ==  mLeftMostArc)  // TODO only for tests
   {
      mLeftMostArc =arc->next;
   }
   if(arc->prev != nullptr)
   {
      arc->prev->next = arc->next;
      arc->next->prev = arc->prev;
   } else {
      arc->next->prev = nullptr;
   }
   // delete from the red black tree
   fixupDelete(arc);
   //delete arc;  TODO not safe if for example arc was allocatd on stack, handle explicitly elsewhere
}

/**
  find arc of the red-black tree directly above the newSite, O(log(n))
  nullptr if no arc found
*/ 
Arc_s* CBeachBST::findArcAbove(Arc_s* newSite, double sweepLineYpos)  
{
   Arc_s* currArc =  mRoot;
   while(currArc != nullptr)
   {
	   Breakpoint_s breakpoints = calculateArcBreakpoints(currArc, sweepLineYpos);
	   if( (newSite->sitePoint).x > breakpoints.rightBreakpoint)
	   {
	      currArc = currArc -> rigthChild;
	   }
	   if( (newSite->sitePoint).x < breakpoints.leftBreakpoint)
	   {
	      currArc = currArc -> leftChild;
	   } else {
		   // we found the arc
	      return currArc;
	   }
   }
   return currArc; 
}

/**
  look for arc in linear time, we will use this before red black tree search is implemented
*/
Arc_s* CBeachBST::findArcAboveInLinearTime(Arc_s* newSite, double sweepLineYpos) // for firsts tests implement simple linear search, later implement red black tree
{
	Arc_s* currArc =  mLeftMostArc;  // look for arc strating with leftmost arc, linear time
	while(currArc != nullptr)
    {
	   Breakpoint_s breakpoints = calculateArcBreakpoints(currArc, sweepLineYpos);
	   if((newSite->sitePoint).x <= breakpoints.rightBreakpoint && (newSite->sitePoint).x >= breakpoints.leftBreakpoint )
	   {
		   // arc found
		   return currArc;
	   }
	   currArc = currArc->next; // else keep looking
    }
	return currArc; // nullptr if no arc found
}

void CBeachBST::fixupInsert(Arc_s* arc) // fixup after inserting a new arc so that all red-black properties of the tree are ensured
{
   //ToDO
}


void CBeachBST::fixupDelete(Arc_s* arc)
{
   //TODO
}

/*
  determine if breakpoints converge based on current breakpoint time evolution vectors
*/
bool CBeachBST::doCurrentBreakPointsConverge(Arc_s* left, Arc_s* middle,Arc_s* right, Point& circleEventPoint , Point& intersectionPoint, double sweelplineY)
{
	 
   if(left== nullptr || middle == nullptr || right == nullptr)
   {
      return false;
   }
   // find the intersection of orthogonal lines starting from middle point left-middle and middle-right -> circle event
   Point A = left->sitePoint;  // get verteces of the triangle of site points
   Point B = middle->sitePoint;
   Point C = right->sitePoint;
   //  calculate breakpoint starting points at the "moment" sweeplineY
   Point AB_Middle = calculateParabolaIntersection (A, B, sweelplineY );
   Point BC_Middle = calculateParabolaIntersection (B, C, sweelplineY );
   
   // solution of the equation:   (l, -t) = A (ABMiddle - BCMiddle)
   // matrix A =  ( ( AB.y, -AB.x ) , ( BC.y, -BC.x ) )
   // adjugate of A =  ( ( -BC.x, AB.x  ), ( -BC.y , AB.y ) )
   // (l, -t) = A^(-1) (leftMiddle - rightMiddle)
   // A_inv = A_adj/ A_det
   // calculate determinant of A to see if A has an inverse
   
    // determine direction of vectors dependent on the development of breakpoints in dependence on sweep line
   auto evolutionVectorpair = calcBrPntEvolutionDirection(middle, sweelplineY); // first element is left breakpoint evoltuion vector, second element is ditto
   Matrix_2x2 A_mat { { {evolutionVectorpair.first.x, evolutionVectorpair.first.y } , { evolutionVectorpair.second.x , evolutionVectorpair.second.y} } };
   double A_det =  calcMatrixDet(A_mat);
   if(approxEquality_n::isEqualZero(A_det)) // colinear breakpoint vectors
   {
      return false;
   }
   Matrix_2x2 A_inv = calcInverseMatrix(A_mat);
   Vector AB_mid_Minus_BC_mid(-AB_Middle.x + BC_Middle.x, -AB_Middle.y + BC_Middle.y);
   Vector vec_l_t = multMatrixVector(A_inv, AB_mid_Minus_BC_mid);  // matrix multiplication step
   double l = vec_l_t.x;  // vector parametrization 
   double t = -vec_l_t.y;
   if(l < 0 || t < 0)
   {
      return false; // intersection is in the past, not in the future
   } 
   //  calculate intersection and circle event point
   Vector DO_mult_t = mult(evolutionVectorpair.first, l);
   intersectionPoint = addVector (AB_Middle, DO_mult_t);
   circleEventPoint.x = intersectionPoint.x;
   circleEventPoint.y = intersectionPoint.y - euclidDist(A, intersectionPoint); // lowest Point of the circle, that includes points A, B, C  
   return true;
}

/*
  calculate how the breakpoint of an arc evolves for a certain sweepline change
  firt element-> left breakpoint vector, second-> right
*/
std::pair<Vector, Vector> CBeachBST::calcBrPntEvolutionDirection(Arc_s* arc, double startSweepLineYpos)
{
  //TODO handle case if next or prev are nullptr
   const double sweepInc1 = 0.1;
   const double sweepInc2 = 1;
   Point right_now = calculateParabolaIntersection (arc->sitePoint, (arc->next)->sitePoint, startSweepLineYpos - sweepInc2 );  // right breakpoints evolution
   Point right_before = calculateParabolaIntersection (arc->sitePoint, (arc->next)->sitePoint, startSweepLineYpos - sweepInc1);
   
   Point left_now = calculateParabolaIntersection ((arc->prev)->sitePoint, arc->sitePoint, startSweepLineYpos - sweepInc2);
   Point left_before = calculateParabolaIntersection ((arc->prev)->sitePoint, arc->sitePoint, startSweepLineYpos - sweepInc1);
   
   return std::make_pair(Vector(left_now.x - left_before.x, left_now.y - left_before.y ), Vector(right_now.x - right_before.x, right_now.y - right_before.y ));
}


/*
  calculate parabola y coord from its focus point p, x coordinate of a poit on parabola and sweepline y
*/
double CBeachBST::calculateParabolaY(Point p, double x , double sweepLine)
{
   double a =  0.5/(  p.y - sweepLine ) ;
   double b =  - p.x/(p.y - sweepLine);
   double c =   (p.x*p.x + p.y*p.y - sweepLine*sweepLine) /  (2* (p.y - sweepLine));
   return  a*x*x + b*x + c;
}

/**
  Point leftSite:  site point of the left parabola
  Point rightSite: site point of the right parabola
  double sweepLineYpos: sweepline y coordinate
  return: x coordinate of the breakpoint of two parabolic arcs
*/
Point CBeachBST::calculateParabolaIntersection (Point leftSite, Point rightSite, double sweepLineYpos )
{
   // calculate two intersection points of the parabolas
   double px= leftSite.x;
   double py = leftSite.y;
   double qx = rightSite.x;
   double qy = rightSite.y;
   double ly = sweepLineYpos;
   
   if(leftSite.y == sweepLineYpos)
   {
      return Point(px, calculateParabolaY(rightSite, px, sweepLineYpos));
   }
   
   if(rightSite.y == sweepLineYpos)
   {
      return Point(qx, calculateParabolaY(leftSite, qx, sweepLineYpos));
   }
   
   //double delta_a = 1 /(2 * (leftSite.y - sweepLineYpos)) - 1/(2 * ( rightSite.y - sweepLineYpos ) ) ;
   double delta_a = 0.5 * ( 1/(leftSite.y - sweepLineYpos) -  1/(rightSite.y - sweepLineYpos)  );
   //double delta_b = -leftSite.x/(leftSite.y - sweepLineYpos) + rightSite.x/(rightSite.y - sweepLineYpos);
   double delta_b =  rightSite.x/(rightSite.y - sweepLineYpos) - leftSite.x/(leftSite.y - sweepLineYpos);
   double delta_c = (1/( 2* (py-ly) * (qy - ly)  )) * ( qy * (px * px + py * py) + ly * ( qx * qx + qy * qy - px * px - py * py) - py * (qy*qy + qx * qx) + ly*ly * (py - qy) ) ; //TODO
 
   // handle case if delta_a is zero
   if( approxEquality_n::isEqualZero(delta_a) )
   {
      double res = -delta_c/delta_b;
      return Point(res, calculateParabolaY(rightSite, res , sweepLineYpos) ) ;
   }
   
   double d = std::sqrt(delta_b * delta_b - 4* delta_a * delta_c);
   double r1 = (- delta_b - d)/(2* delta_a);
   double r2 = (- delta_b + d)/(2* delta_a);
   
   // choose the correct intersection points based on the order of the arcs
   bool is_r1 =  ( ((r1 - px)*(qy - ly)) > ( (r1- qx) * (py - ly) )  );  // criteria that dy/dy if the left arc bigger than that of the right
   double res;
   if(is_r1)
   {
      res = r1;
   } else {
      res = r2;
   }
   return Point(res, calculateParabolaY(rightSite, res , sweepLineYpos) );
}

/**
 return left and right break points of the arc.
 if no left breakpoint exists -> lowest possible double val
 if no right breakpooint->highest possible double val
*/
CBeachBST::Breakpoint_s CBeachBST::calculateArcBreakpoints(Arc_s* arc, double sweepLineYpos)
{
	// highest and lowest possible values by default
   Breakpoint_s result;
   result.leftBreakpoint = -DBL_MAX ; // lowest possible value
   result.rightBreakpoint = DBL_MAX ;
   if(arc->next != nullptr)  // if right neghbor arc exists
   {
      result.rightBreakpoint = calculateParabolaIntersection (arc->sitePoint, (arc->next)->sitePoint, sweepLineYpos ).x;
   }
   if(arc->prev != nullptr)
   {
      result.leftBreakpoint = calculateParabolaIntersection ( (arc->prev)->sitePoint, arc->sitePoint, sweepLineYpos ).x;   
   }
   return result;
}


#if TEST==1

TEST_CASE("insertFirstArcTest") {
    CBeachBST beach;
    Point p(1,1);
    Arc_s arc(p, 1);
    CHECK(beach.isEmpty() == true);
    beach.insertArcAfter(&arc, nullptr);
    CHECK(beach.isEmpty() == false);
}

/*
   before: arc1->righArc
   
   after:  arc1->arc2->righArc
*/
TEST_CASE("insertArcAfterLinkedListTest")
{
   CBeachBST beach;         
	
   Point p0(0,0);
   Arc_s righArc (p0,0);
	
   Point p1(1,1);
   Arc_s arc1(p1, 1);
   arc1.next = &righArc;
   righArc.prev = &arc1;
   
   Point p2(2,2);
   Arc_s arc2(p2, 1);
   
   beach.insertArcAfter(&arc1, nullptr);
   beach.insertArcAfter(&arc2, &arc1);
	
   CHECK(arc1.next == &arc2);
   CHECK(&arc1 == arc2.prev);
   CHECK(arc2.next == &righArc);
   CHECK(&arc2 == righArc.prev); 
}

TEST_CASE("insertArcBeforeNullptrLinkedListTest")
{
   CBeachBST beach;         
		
   Point p1(1,1);
   Arc_s arc1(p1, 1);
   
   Point p2(2,2);
   Arc_s arc2(p2, 1);
   
   beach.insertArcAfter(&arc1, nullptr);
   beach.insertArcAfter(&arc2, &arc1);
	
   CHECK(arc1.next == &arc2);
   CHECK(&arc1 == arc2.prev);
}

/*
  arc1 <-> arc2 <-> arc3
  arc1 <-> arc3
*/
TEST_CASE("deleteArcLinkedListTest")
{
   CBeachBST beach;         
	
   Point p1(1,1);
   Arc_s arc1(p1, 1);
   
   Point p2(2,1);
   Arc_s arc2(p2, 1);
   
   Point p3(3,1);
   Arc_s arc3(p3, 1);
   
   // connect arcs
   beach.insertArcAfter(&arc1, nullptr);
   beach.insertArcAfter(&arc2, &arc1);
   beach.insertArcAfter(&arc3, &arc2);
   
   beach.deleteArc (&arc2);
   CHECK(arc1.next == &arc3);
   CHECK(&arc1 == arc3.prev);
}

/*
   leftArc -> arc1
   leftArc -> arc2 -> arc1
*/

TEST_CASE("insertArcBeforeLinkedListTest")
{
   CBeachBST beach;         
	
   Point p0(0,0);
   Arc_s leftArc (p0,0);
	
   Point p1(1,1);
   Arc_s arc1(p1, 1);
   arc1.prev = &leftArc;
   leftArc.next = &arc1 ;
   
   Point p2(2,2);
   Arc_s arc2(p2, 1); // arc to insert
   
   beach.insertArcBefore(&arc1, nullptr);
   beach.insertArcBefore(&arc2, &arc1);
   
   CHECK(arc2.next == &arc1);
   CHECK(&arc2 == arc1.prev);
   CHECK(leftArc.next == &arc2);
   CHECK(&leftArc == arc2.prev); 
}

TEST_CASE("BreakPointsConvergenceTestNullPtrArc")
{
    CBeachBST beach;    
    Arc_s* nullArc = nullptr;	
    Point p1(1,1);
    Arc_s arc1(p1, 1);   
    Point p2(2,2);
    Arc_s arc2(p2, 1);	
    Point a;
    Point b;	
    bool res = beach.doCurrentBreakPointsConverge(nullArc, &arc1, &arc2, a, b, 1 );
    CHECK( res == false);
}

TEST_CASE("BreakPointsEvolutionConvergenceTestNullNoIntersection")
{
    CBeachBST beach;    
    Point p1(1,0); // points lie on the same line, the orthogonal lines through their midsections do not converge
    Arc_s arc1(p1, 1);
   
    Point p2(2,0);
    Arc_s arc2(p2, 1);
	
    Point p3(3,0);
    Arc_s arc3(p3, 1);
    
    
    arc1.next =  &arc2;  // connect arcs
    arc2.prev =  &arc1;
   
    arc2.next = &arc3;
    arc3.prev = &arc2;
	
    Point a;
    Point b;
	
    bool res = beach.doCurrentBreakPointsConverge(&arc1, &arc2, &arc3, a, b, 0 );
	
    CHECK( res == false);
}


TEST_CASE("BreakPointsEvolutionConvergenceTestIntersectionInTheFuture")
{
    CBeachBST beach;    
    Point p1(1,0); 
    Point p2(2,1);
    Point p3(3,0);
    Arc_s leftMiddlePrev (p1, 0);  // triple of converging arcs
    Arc_s leftMiddle (p2, 1);
    Arc_s leftMiddleNext( p3, 2);
   
    leftMiddlePrev.next =  &leftMiddle;  // connect arcs
    leftMiddle.prev =  &leftMiddlePrev;
   
    leftMiddle.next = &leftMiddleNext;
    leftMiddleNext.prev = &leftMiddle;
	
    Point circleEventPoint;
    Point intersectionPoint;
	
    bool res = beach.doCurrentBreakPointsConverge(&leftMiddlePrev, &leftMiddle, &leftMiddleNext, circleEventPoint, intersectionPoint, 0 );
	
    CHECK(res == true);
    CHECK( approxEquality_n::isEqualZero(intersectionPoint.x - 2) );   // intersection point is (2,0)
    CHECK( approxEquality_n::isEqualZero(intersectionPoint.y) ); 
}

TEST_CASE("ArcBreakpointsIntersectionTestSingleArcWithoutNeighbors")
{
   CBeachBST beach;    
   Point p1(1,0); 
   Arc_s arc1(p1, 1);
   CBeachBST::Breakpoint_s breakpoints = beach.calculateArcBreakpoints(&arc1, 1);
   CHECK( breakpoints.leftBreakpoint == -DBL_MAX);
   CHECK( breakpoints.rightBreakpoint == DBL_MAX);
}

TEST_CASE("ParabolaIntersectionTestBreakpointInTheMiddle")
{
   CBeachBST beach;
   double x1 = 1;
   double x2 = 2;   
   double sweepLineY = -2;  // should be below the y coordinates of both points
   Point p1(x1,0); 
   Point p2(x2,0); 
   double breakpointX = beach.calculateParabolaIntersection (p1, p2, sweepLineY).x;
   //CHECK( beach.isEqualZero(breakpointX - (x1+x2)/2 ) == true );
   CHECK( (x1+x2)/2  == breakpointX);
}

#endif

