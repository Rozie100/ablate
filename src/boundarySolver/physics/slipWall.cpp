#include "slipWall.hpp"
#include <utilities/mathUtilities.hpp>
#include "finiteVolume/compressibleFlowFields.hpp"

void ablate::boundarySolver::physics::SlipWall::Setup(ablate::boundarySolver::BoundarySolver &bSolver) {
    bSolver.RegisterFunction(UpdateBoundaryVel, this, {finiteVolume::CompressibleFlowFields::EULER_FIELD}, {finiteVolume::CompressibleFlowFields::VELOCITY_FIELD});
}

PetscErrorCode ablate::boundarySolver::physics::SlipWall::UpdateBoundaryVel(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg,
                                                                                  const PetscFVCellGeom *boundaryCell, const PetscInt *uOff, PetscScalar *boundaryValues,
                                                                                  const PetscScalar *stencilValues, const PetscInt *aOff, PetscScalar *auxValues, const PetscScalar *stencilAuxValues,
                                                                                  void *ctx) {
    PetscFunctionBeginUser;
    const PetscInt EULER_FIELD = 0;
    const PetscInt VEL = 0;

    PetscReal stencilVel[3];
    PetscReal boundaryVel[dim];

    // map the stencil velocity in normal coord.
    PetscReal stencilDensity = stencilValues[uOff[EULER_FIELD] + finiteVolume::CompressibleFlowFields::RHO];
    for (PetscInt d = 0; d < dim; d++) {
        stencilVel[d] = stencilValues[uOff[EULER_FIELD] + finiteVolume::CompressibleFlowFields::RHOU + d] / stencilDensity;
    }

   
    // update boundary velocities
    PetscScalar boundaryDensity = boundaryValues[uOff[EULER_FIELD] + finiteVolume::CompressibleFlowFields::RHO];
    for (PetscInt d = 0; d < dim; d++) {
        boundaryVel[d] = stencilVel[d];
        boundaryValues[uOff[EULER_FIELD] + finiteVolume::CompressibleFlowFields::RHOU + d] = boundaryVel[d] * boundaryDensity;
        auxValues[aOff[VEL] + d] = boundaryVel[d];
    }

    PetscFunctionReturn(0);
}

#include "registrar.hpp"
REGISTER_WITHOUT_ARGUMENTS(ablate::boundarySolver::BoundaryProcess, ablate::boundarySolver::physics::SlipWall,
                           "updates boundary velocities at the wall using outflow to reconstruct the log law");
