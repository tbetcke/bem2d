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


}
