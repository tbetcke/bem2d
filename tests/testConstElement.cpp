// Test the Point structure
#include "boost/test/unit_test.hpp"
#include <iostream>
#include "../libs/bem2d.h"

BOOST_AUTO_TEST_CASE( constelement )
{

  bem2d::Point p1=bem2d::Point(1,0);
  bem2d::Point p2=bem2d::Point(1,2);

  bem2d::ConstElement e(p1,p2);

    // Functions

  bem2d::Point p3=e.first();
  BOOST_CHECK_CLOSE( p3.x, 1.0, 1E-15);
  BOOST_CHECK_CLOSE( p3.y, 0.0, 1E-15);


  bem2d::Point p4=e.last();
  BOOST_CHECK_CLOSE( p4.x, 1.0, 1E-15);
  BOOST_CHECK_CLOSE( p4.y, 2.0, 1E-15);


  bem2d::Point p5=e.normal(.5);
  BOOST_CHECK_CLOSE( p5.x, 1.0, 1E-15);
  BOOST_CHECK_CLOSE( p5.y, 0.0, 1E-15);

  bem2d::Point p6=e.deriv(.3);
  BOOST_CHECK_CLOSE( p6.x, 0.0, 1E-15);
  BOOST_CHECK_CLOSE( p6.y, 2.0, 1E-15);


}

