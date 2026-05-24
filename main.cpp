#include "dcel.h"
#include "Voronoi_Diagram.h"
#include <vector>
#include <SFML/Graphics.hpp>
#include <vector>

/*
  calculate areas of voronoi cells
*/ 
std::vector<double> voronoi_areas(const std::vector<Point> &p)  // move to main.cpp
{
   CVoronoiDiagram vorDiagram(p); // dependent on the size of p, maybe use heap?
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
   return vorDiagram.calculateAreas();
}

/*
 draw the diagram from three components
*/
void drawVoronoiDiag(std::vector<Point> genPoints, CVoronoiDiagram::BoundaryBox_s bbox, std::vector<Half_edge_s> HalfEdges)
{
	  std::vector<sf::CircleShape> genPointCircles;
	  float pointRadius = 5;
	  for(auto p: genPoints)
	  {
		  sf::CircleShape shape(pointRadius);  // convert gen point coords to circles
		  shape.setOrigin(pointRadius, pointRadius); // set center of cicle relativ to cicle coordinates, needed to center the circle
		  shape.setFillColor(sf::Color::Green);	
		  shape.setPosition( p.x, p.y );	
		  genPointCircles.push_back(shape);
	  }
	  sf::VertexArray lines(sf::LineStrip, 8);
	  int pointCnt = 0;
          for (auto line: bbox.mLineSegments)
          {
            lines[pointCnt++].position = sf::Vector2f(line.mBegin.x, line.mBegin.y);
            lines[pointCnt++].position = sf::Vector2f(line.mEnd.x, line.mEnd.y);
            //lines.color = sf::Color::Green;
          }
	  //TODO  edges 
	  // convert array of edges to vector of line strips with size 2, every edge and its twin -> one line strip
	  std::vector<sf::VertexArray> lineEdgeStrips ;
	  for(auto edge: HalfEdges)
	  {
	    sf::VertexArray lines(sf::LineStrip, 2);
	    lines[0].position = sf::Vector2f( edge.Origin.x, edge.Origin.y  );
	    int twinEdgeInd = edge.TwinEdge;
	    Half_edge_s twinEdge = HalfEdges.at(twinEdgeInd);
	    lines[1].position = sf::Vector2f( twinEdge.Origin.x, twinEdge.Origin.y  );
	    lineEdgeStrips.push_back(lines);
	  }
	  sf::RenderWindow window(sf::VideoMode(800, 600), "Voronoi Tesselation");
          while (window.isOpen()) 
          {
             sf::Event event;
             while (window.pollEvent(event))
             {
               if (event.type ==  sf::Event::Closed)
               {    
                 window.close(); 
               } 
             }
             window.clear();
             for(auto point:  genPointCircles)
             {
               window.draw(point);
             }
             
             for(auto  strip: lineEdgeStrips)  // draw edges intersecting with boundary box
             {
               window.draw(strip);
             }
             window.draw(lines); //draw boundary box
             window.display();            
         }
}


int main(int argc, char * argv[])
{
   // std::vector<Point> vector = { Point(160,220) , Point(220,205), Point(230,225)};
   //std::vector<Point> vector = { Point(160,220), Point(200,300) , Point(220,205),  Point(250, 150), Point(200, 152) };
   std::vector<Point> vector;
   int i = 1;
   if(argc < 2)
   {
      std::cout<<"Command usage is: " << argv[0] << " "<< "PointX1,PointY1 PointX2,PointY2..." << std::endl;
      return 0;
   }
   while (i < argc) 
   {
        //  tokenize chars coordinates and add to vector
        std::string coordinates(argv[i]); 
        auto delim_pos = coordinates.find(",");
        std::string X_coord = coordinates.substr(0, delim_pos);
        std::string Y_coord = coordinates.substr(delim_pos + 1, coordinates.size() - delim_pos - 1 );
        //std::cout<<"X:"<< X_coord <<" " << "Y"<< ": " << Y_coord << std::endl;
        vector.push_back( Point(stod(X_coord), stod(Y_coord) )  );
        i++;
   }
   CVoronoiDiagram vorDiagram(vector);
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
   vorDiagram.intersectEdgesWithBbox();
   drawVoronoiDiag(vorDiagram.getGenPoints(), vorDiagram.getbBox(), vorDiagram.getEdges());
   return 0;
}
