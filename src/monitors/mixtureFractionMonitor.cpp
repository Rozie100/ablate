#include "mixtureFractionMonitor.hpp"

#include <utility>
#include "finiteVolume/compressibleFlowFields.hpp"

ablate::monitors::MixtureFractionMonitor::MixtureFractionMonitor(std::shared_ptr<MixtureFractionCalculator> mixtureFractionCalculator)
    : mixtureFractionCalculator(std::move(std::move(mixtureFractionCalculator))) {}

void ablate::monitors::MixtureFractionMonitor::Register(std::shared_ptr<solver::Solver> solverIn) {
    // Name this monitor
    auto monitorName = solverIn->GetSolverId() + "_mixtureFraction";

    // Define the required fields
    std::vector<std::shared_ptr<domain::FieldDescriptor>> fields{
        std::make_shared<domain::FieldDescription>("zMix", "zMix", domain::FieldDescription::ONECOMPONENT, domain::FieldLocation::SOL, domain::FieldType::FVM),
        std::make_shared<domain::FieldDescription>(ablate::finiteVolume::CompressibleFlowFields::YI_FIELD,
                                                   ablate::finiteVolume::CompressibleFlowFields::YI_FIELD,
                                                   mixtureFractionCalculator->GetEos()->GetSpeciesVariables(),
                                                   domain::FieldLocation::SOL,
                                                   domain::FieldType::FVM),
        std::make_shared<domain::FieldDescription>("densityEnergySource", "energySource", domain::FieldDescription::ONECOMPONENT, domain::FieldLocation::SOL, domain::FieldType::FVM),
        std::make_shared<domain::FieldDescription>(
            "densityYiSource", "densityYiSource", mixtureFractionCalculator->GetEos()->GetSpeciesVariables(), domain::FieldLocation::SOL, domain::FieldType::FVM)};

    // get the required function to compute density
    densityFunction = mixtureFractionCalculator->GetEos()->GetThermodynamicFunction(eos::ThermodynamicProperty::Density, solverIn->GetSubDomain().GetFields());

    // this probe will only work with fV flow with a single mpi rank for now.  It should be replaced with DMInterpolationEvaluate
    auto finiteVolumeSolver = std::dynamic_pointer_cast<ablate::finiteVolume::FiniteVolumeSolver>(solverIn);
    if (!finiteVolumeSolver) {
        throw std::invalid_argument("The MixtureFractionMonitor monitor can only be used with ablate::finiteVolume::FiniteVolumeSolver");
    }
    // get a reference to the tchem reactions instance in the solver
    chemistry = finiteVolumeSolver->FindProcess<ablate::finiteVolume::processes::Chemistry>();

    // call the base function to create the domain
    FieldMonitor::Register(monitorName, solverIn, fields);
}
void ablate::monitors::MixtureFractionMonitor::Save(PetscViewer viewer, PetscInt sequenceNumber, PetscReal time) {
    PetscFunctionBeginUser;
    // get the required fields from the fieldDm and main dm
    const auto& zMixMonitorField = monitorSubDomain->GetField("zMix");
    const auto& yiMonitorField = monitorSubDomain->GetField(ablate::finiteVolume::CompressibleFlowFields::YI_FIELD);
    const auto& energySourceField = monitorSubDomain->GetField("densityEnergySource");
    const auto& densityYiSourceField = monitorSubDomain->GetField("densityYiSource");
    const auto& eulerField = GetSolver()->GetSubDomain().GetField(ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD);
    const auto& densityYiField = GetSolver()->GetSubDomain().GetField(ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD);

    // define a localFVec from the solution dm to compute the source terms
    Vec sourceTermVec = nullptr;
    if (chemistry) {
        DMGetLocalVector(GetSolver()->GetSubDomain().GetDM(), &sourceTermVec) >> utilities::PetscUtilities::checkError;
        VecZeroEntries(sourceTermVec);
        auto fvSolver = std::dynamic_pointer_cast<ablate::finiteVolume::FiniteVolumeSolver>(GetSolver());
        chemistry->AddChemistrySourceToFlow(*fvSolver, sourceTermVec);
    }

    // Get the arrays for the global vectors
    const PetscScalar* solutionFieldArray;
    PetscScalar* monitorFieldArray;
    VecGetArrayRead(GetSolver()->GetSubDomain().GetSolutionVector(), &solutionFieldArray) >> utilities::PetscUtilities::checkError;
    VecGetArray(monitorSubDomain->GetSolutionVector(), &monitorFieldArray) >> utilities::PetscUtilities::checkError;

    // check for the tmpLocalFArray
    const PetscScalar* sourceTermArray = nullptr;
    if (sourceTermVec) {
        VecGetArrayRead(sourceTermVec, &sourceTermArray) >> utilities::PetscUtilities::checkError;
    }

    // March over each cell in the monitorDm
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(monitorSubDomain->GetDM(), 0, &cStart, &cEnd) >> utilities::PetscUtilities::checkError;

    // Get the cells we need to march over
    DMLabel solutionToMonitor;
    DMPlexGetSubpointMap(monitorSubDomain->GetDM(), &solutionToMonitor) >> utilities::PetscUtilities::checkError;

    const PetscInt* monitorToSolution = nullptr;
    IS monitorToSolutionIs = nullptr;
    // if this is a submap, get the monitor to solution
    if (solutionToMonitor) {
        DMPlexGetSubpointIS(monitorSubDomain->GetDM(), &monitorToSolutionIs) >> utilities::PetscUtilities::checkError;
        ISGetIndices(monitorToSolutionIs, &monitorToSolution) >> utilities::PetscUtilities::checkError;
    }

    // save time to get densityFunctionContext
    const auto densityFunctionContext = densityFunction.context.get();

    for (PetscInt monitorPt = cStart; monitorPt < cEnd; ++monitorPt) {
        PetscInt solutionPt = monitorToSolution ? monitorToSolution[monitorPt] : monitorPt;

        // Get the solutionField and monitorField
        const PetscScalar* solutionField = nullptr;
        DMPlexPointGlobalRead(GetSolver()->GetSubDomain().GetDM(), solutionPt, solutionFieldArray, &solutionField) >> utilities::PetscUtilities::checkError;
        PetscScalar* monitorField = nullptr;
        DMPlexPointGlobalRead(monitorSubDomain->GetDM(), monitorPt, monitorFieldArray, &monitorField) >> utilities::PetscUtilities::checkError;

        const PetscScalar* sourceTermField = nullptr;
        if (sourceTermArray) {
            DMPlexPointGlobalRead(GetSolver()->GetSubDomain().GetDM(), solutionPt, sourceTermArray, &sourceTermField) >> utilities::PetscUtilities::checkError;
        }
        // Do not bother in ghost cells
        if (solutionField && monitorField) {
            // compute the density from the solutionPt
            PetscReal density;
            densityFunction.function(solutionField, &density, densityFunctionContext) >> utilities::PetscUtilities::checkError;

            // Copy over and compute yi
            for (PetscInt sp = 0; sp < yiMonitorField.numberComponents; sp++) {
                monitorField[yiMonitorField.offset + sp] = solutionField[densityYiField.offset + sp] / density;
            }

            // Compute mixture fraction
            monitorField[zMixMonitorField.offset] = mixtureFractionCalculator->Calculate(monitorField + yiMonitorField.offset);

            if (sourceTermField) {
                monitorField[energySourceField.offset] = sourceTermField[eulerField.offset + ablate::finiteVolume::CompressibleFlowFields::RHOE];
                PetscArraycpy(monitorField + densityYiSourceField.offset, sourceTermField + densityYiField.offset, densityYiField.numberComponents);
            }
        }
    }

    // cleanup
    if (sourceTermVec) {
        VecRestoreArrayRead(sourceTermVec, &sourceTermArray) >> utilities::PetscUtilities::checkError;
        DMRestoreLocalVector(GetSolver()->GetSubDomain().GetDM(), &sourceTermVec) >> utilities::PetscUtilities::checkError;
    }
    if (monitorToSolutionIs) {
        ISRestoreIndices(monitorToSolutionIs, &monitorToSolution) >> utilities::PetscUtilities::checkError;
    }
    VecRestoreArrayRead(GetSolver()->GetSubDomain().GetSolutionVector(), &solutionFieldArray) >> utilities::PetscUtilities::checkError;
    VecRestoreArray(monitorSubDomain->GetSolutionVector(), &monitorFieldArray) >> utilities::PetscUtilities::checkError;

    // Call the base Save function only after the subdomain global function is updated
    FieldMonitor::Save(viewer, sequenceNumber, time);
    PetscFunctionReturnVoid();
}

#include "registrar.hpp"
REGISTER(ablate::monitors::Monitor, ablate::monitors::MixtureFractionMonitor,
         "This class computes the mixture fraction for each point in the domain and outputs zMix, Yi, and source terms to the hdf5 file",
         ARG(ablate::monitors::MixtureFractionCalculator, "mixtureFractionCalculator", "the calculator used to compute zMix"));