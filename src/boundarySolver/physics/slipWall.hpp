#ifndef ABLATELIBRARY_LODI_INLET_HPP
#define ABLATELIBRARY_LODI_INLET_HPP
#include "boundarySolver/lodi/lodiBoundary.hpp"

#include "boundarySolver/boundaryProcess.hpp"
#include "finiteVolume/processes/navierStokesTransport.hpp"

namespace ablate::boundarySolver::physics {

class SlipWall : public BoundaryProcess {
   private:

    const std::shared_ptr<finiteVolume::processes::PressureGradientScaling> pressureGradientScaling;

   public:
    void Setup(ablate::boundarySolver::BoundarySolver &bSolver) override;

    // update the boundary velocity field using outflow to reconstruct the log-law near the wall (Deardorff, 1970)
    static PetscErrorCode UpdateBoundaryVel(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg, const PetscFVCellGeom *boundaryCell, const PetscInt *uOff,
                                            PetscScalar *boundaryValues, const PetscScalar *stencilValues, const PetscInt *aOff, PetscScalar *auxValues, const PetscScalar *stencilAuxValues,
                                            void *ctx);
};

}  // namespace ablate::boundarySolver::physics
#endif  // ABLATELIBRARY_INLET_HPP
