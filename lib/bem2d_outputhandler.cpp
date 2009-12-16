#include "boost/shared_ptr.hpp"
#include "bem2d_outputhandler.h"
#include "bem2d_outputroutines.h"

#ifdef BEM2DMPI
#include "bem2d_mpi.h"
#endif

namespace bem2d
{

pdvector OutputHandler::RealData(const cvector& vals)
{
  pdvector result(new dvector(vals.size()));
  for (int i=0; i<vals.size(); i++) (*result)[i]=vals[i].real();
  return result;
}

pdvector OutputHandler::ImagData(const cvector& vals)
{
  pdvector result(new dvector(vals.size()));
  for (int i=0; i<vals.size(); i++) (*result)[i]=vals[i].imag();
  return result;
}

pdvector OutputHandler::TurnToRealImag(const cvector& vals)
{
  if (real_)
    {
      return RealData(vals);
    }
  else
    {
      return ImagData(vals);
    }
}


GplotOutput::GplotOutput(int xpts, int ypts, double ax, double bx, double ay, double by, std::string name)
{
  xpts_=xpts;
  ypts_=ypts;
  name_=name;
  set_real(true);
  mesh_=*(MeshGrid(ax,bx,ay,by,xpts_,ypts_));
#ifdef BEM2DMPI
	isroot_=BlacsSystem::Instance()->IsRoot();
#else
	isroot_=true;
#endif
}

void GplotOutput::WriteIncident(const cvector& vals)
{
  if (isroot_) GplotOut(name_+"inc",mesh_,*(TurnToRealImag(vals)),xpts_,ypts_);
}

void GplotOutput::WriteScattered(const cvector& vals)
{
  if (isroot_) GplotOut(name_+"scatt",mesh_,*(TurnToRealImag(vals)),xpts_,ypts_);
}

void GplotOutput::WriteFull(const cvector& vals)
{
  if (isroot_) GplotOut(name_+"full",mesh_,*(TurnToRealImag(vals)),xpts_,ypts_);
}

void GplotOutput::WriteAll(const cvector& valsincident, const cvector& valsscattered, const cvector& valsfull)
{
	if (isroot_){
  GplotOut(name_+"inc",mesh_,*(TurnToRealImag(valsincident)),xpts_,ypts_);
  GplotOut(name_+"scatt",mesh_,*(TurnToRealImag(valsscattered)),xpts_,ypts_);
  GplotOut(name_+"full",mesh_,*(TurnToRealImag(valsfull)),xpts_,ypts_);
	}
}


}

