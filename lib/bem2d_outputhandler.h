#ifndef _OUTPUTHANDLER_H_
#define _OUTPUTHANDLER_H_

#include<vector>
#include "boost/shared_ptr.hpp"
#include <string>
#include "bem2d_defs.h"
#include "bem2d_point.h"

namespace bem2d
{

class OutputHandler
{
public:
        virtual ~OutputHandler() {};
        virtual void WriteIncident(const cvector& vals) {};
        virtual void WriteScattered(const cvector& vals) {};
        virtual void WriteFull(const cvector& vals) {};
        virtual void WriteAll(const cvector& valsincident, const cvector& valsscattered, const cvector& valsfull) {};
        inline const std::vector<Point>& mesh() const {
                return mesh_;
        }
        inline void set_real(bool real) {
                real_=real;
        }

protected:
        pdvector RealData(const cvector& vals);
        pdvector ImagData(const cvector& vals);
        pdvector TurnToRealImag(const cvector& vals);
        std::vector<Point> mesh_;
        bool real_;
        bool isroot_;
};

typedef boost::shared_ptr<OutputHandler> pOutputHandler;

class GplotOutput: public OutputHandler
{
public:
        GplotOutput(int xpts, int ypts, double ax, double bx, double ay, double by, std::string name);
        void WriteIncident(const cvector& vals);
        void WriteScattered(const cvector& vals);
        void WriteFull(const cvector& vals);
        void WriteAll(const cvector& valsincident, const cvector& valsscattered, const cvector& valsfull);
private:
        std::string name_;
        int xpts_;
        int ypts_;
};

}

#endif // _OUTPUTHANDLER_H_

