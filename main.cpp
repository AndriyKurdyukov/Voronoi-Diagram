#include "dcel.h"
#include "Voronoi_Diagram.h"
#include <vector>
#include <SFML/Graphics.hpp>
#include <vector>
#include <fstream>

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
void drawVoronoiDiag(std::vector<Point> genPoints, CVoronoiDiagram::BoundaryBox_s bbox, std::vector<Half_edge_s> HalfEdges, bool isDrawAreas, std::vector<double> cellAreas)
{
   sf::Font font;
   if (!font.loadFromFile("FreeMonoBold.ttf")) {
      return;
   }
   std::string scoreString =  "";          // 3. Set up the sf::Text object
   sf::Text text;
   text.setFont(font);
   text.setString(scoreString);
   text.setCharacterSize(15);            // Size in pixels
   text.setFillColor(sf::Color::Red);  // Text color
   text.setPosition(350.f, 350.f);      
          
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
      int areaIndex = 0;
      for(auto point:  genPointCircles)
      {
         window.draw(point);
         if(isDrawAreas)
         {
            auto pos = point.getPosition();
            text.setString( std::to_string( static_cast<int>(cellAreas.at(areaIndex))) );
            text.setPosition(pos);  
            window.draw(text);   
            areaIndex++;
         }        
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
   //std::vector<Point> vector = { Point(160,220), Point(200,300) , Point(220,205),  Point(250, 150), Point(200, 152) };
   std::vector<Point> vector;  // vector of generator points
   int i = 1;
   if(argc < 2)
   {
      std::cout<<"Command usage is: \n"  <<argv[0] << " [-a] "<< "filename" << "\n" <<
                                         argv[0]   << " -v " << std::endl;
      return 0;
   }
   if(std::string(argv[1]) == "-a") // calculate areas
   {
      i = 2;
   }
   std::string version = "untagged_version";
   if(std::string(argv[1]) == "-v") // output version
   {
     std::cout<< version << std::endl;
     return 0;
   }
   std::ifstream myfile;
   std::string line;
   myfile.open(argv[i]);
   if ( myfile.is_open() ) {
      while ( std::getline (myfile, line) ) 
      {
         auto delim_pos = line.find(",");
         assert(delim_pos != std::string::npos && "coordinates should be comma delimited");
         std::string X_coord = line.substr(0, delim_pos);
         std::string Y_coord = line.substr(delim_pos + 1, line.size() - delim_pos - 1 );
         vector.push_back( Point(stod(X_coord), stod(Y_coord) )  );
      }
   } else {
     std::cout <<"File not found!"<< std::endl;
     return 0;
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
   std::vector<double> cellAreas;
   bool isDrawArea = false;
   if(std::string(argv[1]) == "-a") // calculate areas
   {
      isDrawArea = true;
      cellAreas = vorDiagram.calculateAreas();
   }
   vorDiagram.constructBbox();
   vorDiagram.intersectEdgesWithBbox();
   drawVoronoiDiag(vorDiagram.getGenPoints(), vorDiagram.getbBox(), vorDiagram.getEdges(), isDrawArea, cellAreas);
   return 0;
}
