// Test the Point structure
#include "boost/test/unit_test.hpp"
#include "../libs/bem2d.h"

BOOST_AUTO_TEST_CASE( point )
{
  bem2d::Point p1=bem2d::Point(1.3,2.5);
  bem2d::Point p2=bem2d::Point(3,4);

  // Stored data

  BOOST_CHECK_EQUAL( p1.x, 1.3 );
  BOOST_CHECK_EQUAL( p1.y, 2.5);
  
  // Functions

  BOOST_CHECK_CLOSE( length(p2), 5.0, 1E-15);

  bem2d::Point p3=bem2d::normalize(p2);
  BOOST_CHECK_CLOSE( p3.x, 3./5, 1E-15);
  BOOST_CHECK_CLOSE( p3.y, 4./5, 1E-15);

  // Operators

  bem2d::Point p4=p1+p2;
  BOOST_CHECK_CLOSE( p4.x, p1.x+p2.x,1E-15);
  BOOST_CHECK_CLOSE( p4.y, p1.y+p2.y,1E-15);

  bem2d::Point p5=p1-p2;
  BOOST_CHECK_CLOSE( p5.x, p1.x-p2.x,1E-15);
  BOOST_CHECK_CLOSE( p5.y, p1.y-p2.y,1E-15);

  bem2d::Point p6=3.5*p1;
  BOOST_CHECK_CLOSE( p6.x, 3.5*p1.x,1E-15);
  BOOST_CHECK_CLOSE( p6.y, 3.5*p1.y,1E-15);

  bem2d::Point p7=p1*3.5;
  BOOST_CHECK_CLOSE( p7.x, p6.x, 1E-15);
  BOOST_CHECK_CLOSE( p7.y, p6.y, 1E-15);

}
